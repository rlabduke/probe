import libtbx.load_env

Import("env_base", "env_etc")

env = env_base.Clone(
  LIBS=env_etc.libm)
if (libtbx.manual_date_stamp < 20090819):
  # XXX backward compatibility 2009-08-19
  env.Replace(CCFLAGS=env_etc.ccflags_base)
if (env_etc.compiler != "win32_cl"):
  env.Replace(LINK=env_base["CC"])

exe = env.Program(
  target=["#probe/exe/probe"],
  source=[
    "hybrid_36_c.c",
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
