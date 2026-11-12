"""Execution aliases preserve the canonical fixed-field identity and geometry."""
import re

FIXED_ID = re.compile(r'TJL[45]_\d{3,}')


def controller_name(field_id):
    if not FIXED_ID.fullmatch(field_id):
        raise ValueError('Expected a canonical fixed Trojan field ID')
    return 'MP_'+field_id


def canonical_name(execution_name):
    field_id = execution_name[3:] if execution_name.startswith('MP_') else execution_name
    if not FIXED_ID.fullmatch(field_id):
        raise ValueError('Not a fixed Trojan field name')
    return field_id


def is_fixed_name(name):
    return bool(FIXED_ID.fullmatch(name[3:] if name.startswith('MP_') else name))
