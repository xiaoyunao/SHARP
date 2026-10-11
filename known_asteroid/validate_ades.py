#!/usr/bin/env python3
"""Validate one ADES 2017 PSV block locally with the pinned submission XSD."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

SCHEMA = Path(__file__).resolve().parents[1] / 'heliolincrr' / 'ades_submit.xsd'


def validate_text(text: str, schema_path: Path = SCHEMA) -> dict:
    from lxml import etree

    lines = [line for line in text.splitlines() if line.strip()]
    if not lines or lines[0].replace(' ', '') != '#version=2017':
        raise ValueError('Expected ADES 2017 header')
    root = etree.Element('ades', version='2017')
    block = etree.SubElement(root, 'obsBlock')
    context = etree.SubElement(block, 'obsContext')
    data = etree.SubElement(block, 'obsData')
    parent = None
    columns = None
    nrows = 0
    for line in lines[1:]:
        if line.startswith(('#', '!')):
            if columns is not None:
                raise ValueError('Only one ADES block is supported')
            parts = line[1:].strip().split(None, 1)
            if not parts:
                raise ValueError('Empty ADES metadata')
            name, value = parts[0], parts[1] if len(parts) == 2 else ''
            if line.startswith('#'):
                parent = etree.SubElement(context, name)
                if value:
                    parent.text = value
            else:
                if parent is None:
                    raise ValueError('ADES metadata child without parent')
                etree.SubElement(parent, name).text = value
            continue
        values = [v.strip() for v in line.split('|')]
        if columns is None:
            columns = values
            if not all(columns) or len(set(columns)) != len(columns):
                raise ValueError('Empty or duplicate PSV columns')
            continue
        if len(values) != len(columns):
            raise ValueError(f'Wrong PSV column count at observation {nrows + 1}')
        optical = etree.SubElement(data, 'optical')
        for name, value in zip(columns, values):
            if value:
                etree.SubElement(optical, name).text = value
        nrows += 1
    if not nrows:
        raise ValueError('No ADES observations')
    schema = etree.XMLSchema(etree.parse(str(schema_path)))
    if not schema.validate(root):
        raise ValueError('ADES schema validation failed:\n' + '\n'.join(f'{e} ({e.path})' for e in schema.error_log))
    return {'status': 'local_validated', 'observations': nrows,
            'schema_sha256': hashlib.sha256(schema_path.read_bytes()).hexdigest(),
            'psv_sha256': hashlib.sha256(text.encode()).hexdigest(), 'network_requests': 0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--psv', type=Path)
    source.add_argument('--stdin', action='store_true')
    parser.add_argument('--schema', type=Path, default=SCHEMA)
    args = parser.parse_args()
    text = sys.stdin.read() if args.stdin else args.psv.read_text()
    print(json.dumps(validate_text(text, args.schema), sort_keys=True))


if __name__ == '__main__':
    main()
