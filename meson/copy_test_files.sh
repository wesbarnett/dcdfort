#!/bin/sh

if [[ -v MESON_SOURCE_ROOT ]]; then
    cp "${MESON_SOURCE_ROOT}/src/tests/test.ndx" "${MESON_BUILD_ROOT}/"
    cp "${MESON_SOURCE_ROOT}/src/tests/test.dcd" "${MESON_BUILD_ROOT}/"
fi
