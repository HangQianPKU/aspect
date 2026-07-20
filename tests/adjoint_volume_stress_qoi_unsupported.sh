#!/bin/bash

set -e

output_file="$(mktemp)"
cat > "${output_file}"

grep -q "weak-zone mask and normal/tangent direction interface" "${output_file}"
echo "unsupported plate-boundary stress QOI failed with the expected message"
