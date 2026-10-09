#!/usr/bin/env python3
"""
strip_metadata.py -- remove embedded metadata from every file in a directory.

Usage
    strip_metadata.py DIR                  # cleaned copies go to DIR_stripped/
    strip_metadata.py DIR -o OUTDIR        # cleaned copies go to OUTDIR/
    strip_metadata.py DIR --in-place       # modify the files in DIR itself
    strip_metadata.py DIR --dry-run        # only report what would be done

Options
    --no-recursive       only the top level of DIR (default: all subfolders)
    --keep-orientation   keep the EXIF Orientation tag so photos are not
                         displayed rotated (everything else is removed)
    --keep-xattrs        keep macOS extended attributes (default: remove)

What is removed
  * Images, video, audio, PDF, EPS and every other format ExifTool can
    write (asked from `exiftool -listwf` at run time, so .jpg/.JPG/.jpeg,
    .heic, .png, .tif, .mov, .mp4, camera raw, ... are all covered, in any
    letter case): all EXIF, GPS, XMP, IPTC, ICC, maker notes, comments,
    thumbnails, etc. (`exiftool -all=`).
  * PDF: ExifTool's PDF edits are reversible (the old metadata stays in the
    file as an earlier revision), so if `qpdf` is installed the file is then
    rewritten without old revisions. Without qpdf a warning is printed.
  * Microsoft Office (.docx .xlsx .pptx and macro/template variants): the
    document properties (author, last modified by, company, manager,
    created/modified times, title, subject, keywords, template, editing time,
    custom properties) are cleared, and author names/initials on comments
    and tracked changes are blanked.
  * OpenDocument (.odt .ods .odp .odg ...): meta.xml is cleared.
  * macOS extended attributes on every file (e.g. "Where from" download
    URLs, quarantine flags, Finder tags), unless --keep-xattrs.
  * Finder clutter (.DS_Store, Icon\\r) is not copied (or, with --in-place,
    is deleted).

What is not changed (reported at the end)
  * Formats with no handler, e.g. netCDF/HDF attributes, .mat files, legacy
    .doc/.xls/.ppt, zip archives. Plain-text files (.txt .csv .m .py ...)
    have no embedded metadata block and are copied unchanged.
  * File-system dates (creation/modification times) are not altered.
  * Text a person typed into a document (names in the body, etc.).

Requires ExifTool (https://exiftool.org); qpdf is recommended for PDFs.
Written 9 October 2026 by Ben Kravitz (bkravitz@iu.edu or
ben.kravitz.work@gmail.com) with Claude; replaces test.py, which skipped
.jpg files and upper-case extensions.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

OOXML = {'.docx', '.docm', '.dotx', '.dotm', '.xlsx', '.xlsm', '.xltx', '.xltm',
         '.pptx', '.pptm', '.potx', '.potm', '.ppsx', '.ppsm'}
ODF = {'.odt', '.ods', '.odp', '.odg', '.odf', '.ott', '.ots', '.otp', '.otg'}
TEXT = {'.txt', '.csv', '.tsv', '.md', '.m', '.py', '.r', '.f', '.f90', '.c', '.h',
        '.cpp', '.java', '.js', '.sh', '.tex', '.bib', '.json', '.yaml', '.yml',
        '.ncl', '.ini', '.cfg', '.log', '.dat'}
FINDER_JUNK = {'.DS_Store', 'Icon\r'}
# ExifTool groups that describe the file itself rather than embedded metadata
NEUTRAL_GROUPS = {'ExifTool', 'System', 'File', 'Composite', 'JFIF', 'PNG',
                  'PNG-IHDR', 'QuickTime', 'PDF', 'GIF', 'RIFF', 'MPEG', 'BMP', 'ZIP'}
# TIFF image-structure tags (needed to decode the image; not personal metadata)
TIFF_STRUCTURE = {'ImageWidth', 'ImageHeight', 'BitsPerSample', 'Compression',
                  'PhotometricInterpretation', 'StripOffsets', 'SamplesPerPixel',
                  'RowsPerStrip', 'StripByteCounts', 'XResolution', 'YResolution',
                  'PlanarConfiguration', 'ResolutionUnit', 'TileWidth', 'TileLength',
                  'TileOffsets', 'TileByteCounts', 'SubfileType', 'Predictor',
                  'ExtraSamples', 'SampleFormat', 'ColorMap', 'FillOrder'}

EMPTY_CORE = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
              '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/'
              'package/2006/metadata/core-properties" '
              'xmlns:dc="http://purl.org/dc/elements/1.1/" '
              'xmlns:dcterms="http://purl.org/dc/terms/" '
              'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
              'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"/>')
EMPTY_CUSTOM = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/'
                '2006/custom-properties" xmlns:vt="http://schemas.openxmlformats.org/'
                'officeDocument/2006/docPropsVTypes"/>')
EMPTY_ODF_META = ('<?xml version="1.0" encoding="UTF-8"?>\n'
                  '<office:document-meta xmlns:office="urn:oasis:names:tc:opendocument:'
                  'xmlns:office:1.0" office:version="1.2"><office:meta/>'
                  '</office:document-meta>')
APP_FIELDS = ('Company', 'Manager', 'HyperlinkBase', 'Template', 'TotalTime',
              'Application', 'AppVersion')
RAW = {'.arw', '.cr2', '.cr3', '.nef', '.nrw', '.orf', '.raf', '.rw2', '.pef', '.srw',
       '.erf', '.3fr', '.iiq', '.mos', '.mef', '.x3f'}
# names of people attached to comments and tracked changes; only these attributes
# are touched (never generic ones like name=, which also label fonts and styles)
PERSON_ATTRS = re.compile(rb'(\s(?:[A-Za-z0-9]+:)?(?:author|initials|userId|providerId))="[^"]*"')
AUTHOR_ELEM = re.compile(rb'<author>[^<]*</author>')
REVISION_DATE = re.compile(rb'\sw:date="[^"]*"')
# parts that list people by name (PowerPoint comment authors, Excel threaded-comment persons)
PEOPLE_PARTS = re.compile(r'^(ppt/commentAuthors\.xml|ppt/authors\.xml|xl/persons/[^/]+\.xml)$')
PEOPLE_NAMES = re.compile(rb'(\s(?:name|displayName))="[^"]*"')


def run(cmd):
    """Run a command (list form: no shell), returning (ok, combined output)."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode == 0, (p.stdout + p.stderr).strip()


def exiftool_writable():
    ok, out = run(['exiftool', '-listwf'])
    if not ok:
        sys.exit('ExifTool not found. Install it (e.g. "sudo port install p5-image-exiftool" '
                 'or from https://exiftool.org) and try again.')
    words = out.split(':', 1)[-1].split()
    return {'.' + w.lower() for w in words}


def strip_exiftool(path, keep_orientation):
    cmd = ['exiftool', '-q', '-q', '-overwrite_original', '-all=']
    if keep_orientation:
        cmd += ['-tagsfromfile', '@', '-Orientation']
    ok, out = run(cmd + [str(path)])  # absolute path, so never read as an option
    if not ok:
        return ok, out
    # In TIFF-based files (TIFF, DNG, camera raw) descriptive tags such as Artist,
    # Make, Software and dates share IFD0 with the image structure, and ExifTool
    # will not delete IFD0 as a group; delete each remaining one by name.
    extra = [t for t in remaining_tags(path)
             if t.split(':', 1)[0] in ('IFD0', 'IFD1', 'SubIFD')
             and not (keep_orientation and t.endswith(':Orientation'))]
    if extra:
        ok, out = run(['exiftool', '-q', '-q', '-overwrite_original']
                      + [f'-{t}=' for t in extra] + [str(path)])
    return ok, out


def clean_pdf(path, have_qpdf):
    if not have_qpdf:
        return 'qpdf not installed: earlier PDF revisions with the old metadata remain recoverable'
    tmp = path.with_name(path.name + '.qpdftmp')
    ok, out = run(['qpdf', '--linearize', str(path), str(tmp)])
    if ok or tmp.exists():
        tmp.replace(path)
        return None
    return f'qpdf failed ({out}); old PDF revisions may remain'


def rewrite_zip(path, transform):
    """Rewrite a zip-based document, passing each member through transform(name, bytes).
    Member order and compression are preserved (ODF needs 'mimetype' first, stored)."""
    tmp = path.with_name(path.name + '.ziptmp')
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, 'w') as zout:
        for info in zin.infolist():
            data = transform(info.filename, zin.read(info.filename))
            if data is None:
                continue
            new = zipfile.ZipInfo(info.filename, date_time=(1980, 1, 1, 0, 0, 0))
            new.compress_type = info.compress_type
            new.external_attr = info.external_attr
            zout.writestr(new, data)
    tmp.replace(path)


def ooxml_transform(name, data):
    if name == 'docProps/core.xml':
        return EMPTY_CORE.encode()
    if name == 'docProps/custom.xml':
        return EMPTY_CUSTOM.encode()
    if name == 'docProps/app.xml':
        for field in APP_FIELDS:
            data = re.sub(rb'<%s>[^<]*</%s>' % (field.encode(), field.encode()),
                          b'<%s></%s>' % (field.encode(), field.encode()), data)
        return data
    if name.endswith('.xml') and not name.startswith('docProps/'):
        # names/initials on comments and tracked changes
        data = PERSON_ATTRS.sub(rb'\1=""', data)
        data = REVISION_DATE.sub(b'', data)  # optional w:date on comments/revisions
        data = AUTHOR_ELEM.sub(b'<author></author>', data)
        if PEOPLE_PARTS.match(name):
            data = PEOPLE_NAMES.sub(rb'\1=""', data)
    return data


def odf_transform(name, data):
    if name == 'meta.xml':
        return EMPTY_ODF_META.encode()
    return data


def remaining_tags(path):
    """Embedded-metadata tags ExifTool still finds (excluding file-structure groups)."""
    ok, out = run(['exiftool', '-j', '-G1', '-a', str(path)])
    if not ok and not out.startswith('['):
        return []
    try:
        tags = json.loads(out)[0]
    except (ValueError, IndexError):
        return []
    left = []
    for k in tags:
        if ':' not in k:
            continue
        group, tag = k.split(':', 1)
        if group in NEUTRAL_GROUPS or group.startswith(('Track', 'PNG', 'ZIP')):
            continue
        if group in ('IFD0', 'IFD1', 'SubIFD') and tag in TIFF_STRUCTURE:
            continue
        left.append(k)
    return sorted(left)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir', type=Path)
    where = ap.add_mutually_exclusive_group()
    where.add_argument('-o', '--output', type=Path)
    where.add_argument('--in-place', action='store_true')
    ap.add_argument('--no-recursive', action='store_true')
    ap.add_argument('--keep-orientation', action='store_true')
    ap.add_argument('--keep-xattrs', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    src = a.dir.expanduser().resolve()
    if not src.is_dir():
        sys.exit(f'{src} is not a directory')
    writable = exiftool_writable()
    have_qpdf = shutil.which('qpdf') is not None

    if a.in_place:
        dst = src
    else:
        dst = (a.output.expanduser().resolve() if a.output
               else src.with_name(src.name + '_stripped'))
        if dst == src or src in dst.parents:
            sys.exit('The output directory must not be inside the input directory.')
        if dst.exists() and any(dst.iterdir()) and not a.dry_run:
            sys.exit(f'{dst} already exists and is not empty; choose another -o.')

    files = sorted(p for p in (src.iterdir() if a.no_recursive else src.rglob('*'))
                   if p.is_file() and not p.is_symlink())
    report = {'stripped': [], 'office': [], 'odf': [], 'text': [], 'unsupported': [],
              'junk': [], 'errors': [], 'warnings': [], 'leftover': []}

    for f in files:
        rel = f.relative_to(src)
        if f.name in FINDER_JUNK:
            report['junk'].append(rel)
            if a.in_place and not a.dry_run:
                f.unlink()
            continue
        ext = f.suffix.lower()
        kind = ('office' if ext in OOXML else 'odf' if ext in ODF else
                'stripped' if ext in writable else 'text' if ext in TEXT else 'unsupported')
        report[kind].append(rel)
        if a.dry_run:
            continue
        target = dst / rel
        if not a.in_place:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, target)  # contents only: no xattrs, no Finder info
        try:
            if kind == 'stripped':
                ok, out = strip_exiftool(target, a.keep_orientation)
                if not ok:
                    report['errors'].append(f'{rel}: ExifTool failed: {out}')
                    if not a.in_place:
                        target.unlink()  # never leave an unstripped copy behind
                    continue
                if ext == '.pdf':
                    warn = clean_pdf(target, have_qpdf)
                    if warn:
                        report['warnings'].append(f'{rel}: {warn}')
            elif kind == 'office':
                rewrite_zip(target, ooxml_transform)
            elif kind == 'odf':
                rewrite_zip(target, odf_transform)
            if kind in ('stripped', 'office', 'odf'):
                left = remaining_tags(target)
                if left:
                    report['leftover'].append(f'{rel}: {", ".join(left[:8])}'
                                              + (' ...' if len(left) > 8 else ''))
        except (zipfile.BadZipFile, OSError) as e:
            report['errors'].append(f'{rel}: {e}')
            if not a.in_place and target.exists():
                target.unlink()
            continue
        if not a.keep_xattrs:
            run(['xattr', '-c', str(target)])

    # ---- summary
    print(f'{"Would process" if a.dry_run else "Processed"} {len(files)} files in {src}')
    if not a.in_place and not a.dry_run:
        print(f'Cleaned copies written to {dst}')
    labels = [('stripped', 'metadata stripped (ExifTool)'),
              ('office', 'Office documents cleaned'), ('odf', 'OpenDocument files cleaned'),
              ('text', 'plain text (no embedded metadata)'),
              ('unsupported', 'NOT CHANGED: no handler for this type'),
              ('junk', 'Finder files (.DS_Store/Icon) skipped')]
    for key, label in labels:
        if report[key]:
            print(f'  {len(report[key]):6d}  {label}')
    if report['unsupported']:
        exts = sorted({p.suffix.lower() or '(no extension)' for p in report['unsupported']})
        print('          unsupported types: ' + ' '.join(exts))
    for key, title in (('warnings', 'Warnings'), ('leftover', 'Tags still present (check)'),
                       ('errors', 'ERRORS (file not cleaned' +
                        ('' if a.in_place else '; no copy written') + ')')):
        if report[key]:
            print(f'{title}:')
            for line in report[key]:
                print('  ' + line)
    return 1 if report['errors'] else 0


if __name__ == '__main__':
    sys.exit(main())
