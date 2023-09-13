#!/bin/bash
case "$PWD" in
  *\ * )
    echo "Cannot install into $PWD. Directory must not contain spaces."
    exit 1
    ;;
esac

source bin/activate
conda-unpack
