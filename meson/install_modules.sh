#!/bin/sh

# Installs Fortran module files into include directory

if [ -n "${MESON_INSTALL_DESTDIR_PREFIX}" ]; then
    install -v -Dm644 "${MESON_BUILD_ROOT}"/dcdfort@sha/*.mod -t "${MESON_INSTALL_DESTDIR_PREFIX}"/include/
fi
