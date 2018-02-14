#!/bin/sh

# Installs Fortran module files into include directory
project="dcdfort"

if [[ -v MESON_INSTALL_DESTDIR_PREFIX ]]; then
    install -Dm644 "${MESON_BUILD_ROOT}"/"$project"@sha/*.mod -t "${MESON_INSTALL_DESTDIR_PREFIX}"/include/
fi
