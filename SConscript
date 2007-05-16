import libtbx.load_env

Import("env_base", "env_etc")

env = env_base.Copy(
  CCFLAGS=env_etc.ccflags_base,
  LIBS=env_etc.libm
)

exe = env.Program(
  target=["#probe/exe/probe"],
  source=[
    "probe.c",
    "dots.c",
    "abin.c",
    "readPDBrecs.c",
    "geom3d.c",
    "utility.c",
    "select.c",
    "parse.c",
    "atomprops.c",
    "stdconntable.c",
    "autobondrot.c"])
