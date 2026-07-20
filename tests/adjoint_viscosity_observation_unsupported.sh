#!/bin/bash

set -e

output_file="$(mktemp)"
cat > "${output_file}"

grep -q "direct viscosity observation objectives are parsed" "${output_file}"
echo "unsupported direct viscosity observation failed with the expected message"
