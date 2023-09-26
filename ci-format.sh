#!/usr/bin/env bash

set -eu

type -f clang-format-12 || {
  echo "install clang-format-12"
  exit 1
}

fmt="clang-format-12 -i {} && echo {}"
set -x
for folder in src ibst; do
find \
  $folder/ -name '*.cxx' -exec sh -c "$fmt" \; \
  -or -name '*.hpp' -exec sh -c "$fmt" \; \
  -or -name '*.h' -exec sh -c "$fmt" \;
done
