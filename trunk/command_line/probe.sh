#! /bin/sh
# LIBTBX_SET_DISPATCHER_NAME phenix.probe

if [ -n "$LIBTBX__VALGRIND_FLAG__" ]; then
  exec $LIBTBX_VALGRIND $LIBTBX_BUILD/probe/exe/probe "$@"
elif [ $# -eq 0 ]; then
  exec $LIBTBX_BUILD/probe/exe/probe
else
  exec $LIBTBX_BUILD/probe/exe/probe "$@"
fi
