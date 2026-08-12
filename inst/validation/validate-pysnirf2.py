#!/usr/bin/env python3
"""Validate SNIRF files with pysnirf2 0.7.3 under NumPy >= 2."""

import json
import sys

import numpy as np

if not hasattr(np, "string_"):
    np.string_ = np.bytes_

from pysnirf2 import validateSnirf


def issue_dict(issue):
    return {
        "location": str(issue.location),
        "severity": int(issue.severity),
        "name": str(issue.name),
        "message": str(issue.message),
    }


output = []
for filename in sys.argv[1:]:
    result = validateSnirf(filename)
    output.append(
        {
            "file": filename,
            "valid": bool(result),
            "errors": [issue_dict(x) for x in result.errors],
            "warnings": [issue_dict(x) for x in result.warnings],
        }
    )

print(json.dumps(output, sort_keys=True, separators=(",", ":")))
