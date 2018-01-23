#!/bin/sh

# Installs Fortran module files into include directory

install -v -Dm644 ${MESON_BUILD_ROOT}/dcdfort@sha/*.mod -t ${MESON_INSTALL_PREFIX}/include/
