
# ---------------------------------------------------------------------------
# marley-config: static body
#
# Everything above this line is generated at build time as a block of
# KEY='value' assignments; everything below ships verbatim. The two halves
# are concatenated, not substituted. See design spec 6.2.
#
# The division carries the whole design: MARLEY-owned paths are computed
# HERE, from this script's own location, because the build system that wrote
# the header cannot know where the user will install. External dependency
# flags are baked into the HEADER, because this script cannot know what the
# build actually found. Each half is incapable of making the other's
# mistake, which is why a baked build-tree -L or source-tree -I can no
# longer appear in an installed marley-config.
# ---------------------------------------------------------------------------

# Fail loudly, in one place, if the generated header is missing or was only
# partly written. This matters because an unsubstituted placeholder is not a
# shell syntax error: without this check the script would run to completion
# and print quietly wrong flags.
: "${version:?marley-config: generated header missing or incomplete}"
: "${cxx_std:?marley-config: generated header missing or incomplete}"
: "${use_root:?marley-config: generated header missing or incomplete}"
: "${hepmc3_mode:?marley-config: generated header missing or incomplete}"
: "${gsl_mode:?marley-config: generated header missing or incomplete}"
: "${libdir_rel:?marley-config: generated header missing or incomplete}"
: "${bindir_rel:?marley-config: generated header missing or incomplete}"
: "${includedir_rel:?marley-config: generated header missing or incomplete}"

# Self-locate. BASH_SOURCE rather than $0, matching the other shipped
# scripts: the two differ when a script is sourced rather than executed.
here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
prefix="$( cd "${here}/.." && pwd )"

bindir="${prefix}/${bindir_rel}"
libdir="${prefix}/${libdir_rel}"
includedir="${prefix}/${includedir_rel}"

# --- compiler flags --------------------------------------------------------

cflags="-I${includedir} -std=${cxx_std}"

# The vendored HepMC3 headers are nested under the MARLEY include root, so
# the -I that reaches them is derived here like any other MARLEY-owned path.
# An external HepMC3 lives somewhere only the build system knew about.
if [ "$hepmc3_mode" = "builtin" ]; then
  cflags="${cflags} -I${includedir}/marley/builtin"
else
  cflags="${cflags} ${ext_hepmc3_cflags}"
fi

# -DUSE_ROOT is not optional on a ROOT-enabled build: it gates a data member
# in EventFileReader.hh, so a consumer compiling without it disagrees with
# the installed library about the layout of a public class.
if [ "$use_root" = "yes" ]; then
  cflags="${cflags} ${root_cflags} -DUSE_ROOT"
fi

# GSL is deliberately absent from both cflags and libs. It is a private
# dependency -- used only by src/marley_gsl.cc, never named in an installed
# header -- so a consumer has no reason to compile or link against it. The
# gsl_* keys remain in the generated header for diagnosis.

# --- linker flags ----------------------------------------------------------

libs="-L${libdir} -lMARLEY"

if [ "$hepmc3_mode" = "builtin" ]; then
  libs="${libs} -lMARLEYHepMC3"
else
  libs="${libs} ${ext_hepmc3_libs}"
fi

if [ "$use_root" = "yes" ]; then
  libs="${libs} ${root_libs}"
fi

# The rpath is carried here but NOT in marley.pc. The audiences differ: a
# hand-written Makefile has no other way to make its binary runnable, while
# a packager's toolchain manages rpath itself and an unexpected one is an
# obstacle.
libs="${libs} ${new_dtags} -Wl,-rpath,${libdir}"

# --- command line ----------------------------------------------------------

usage_msg="Usage: marley-config [--cflags] [--libs] [--cxx-std] [--prefix]"
usage_msg="${usage_msg} [--bindir] [--libdir] [--includedir] [--use-root]"
usage_msg="${usage_msg} [--version] [--git-revision] [--help]"

if [ $# -eq 0 ]; then
  echo "${usage_msg}"
  exit 1
fi

for arg in "$@"; do
  if [ "$arg" = "--help" ]; then
    echo "Usage: marley-config [options]"
    echo ""
    echo "  --cflags          Print compiler flags and header path"
    echo "  --libs            Print linker flags for the MARLEY libraries"
    echo "  --cxx-std         Print the C++ standard MARLEY was built with"
    echo "  --prefix          Print the installation prefix"
    echo "  --bindir          Print the executable directory"
    echo "  --libdir          Print the library directory"
    echo "  --includedir      Print the header directory"
    echo "  --use-root        Print whether ROOT support is enabled"
    echo "  --version         Print the MARLEY version"
    echo "  --git-revision    Print the MARLEY git revision"
    echo "  --help            Print this message"
    exit 2
  fi
done

out=""

for arg in "$@"; do
  case $arg in
    --cflags)
      out="$out $cflags"
    ;;
    --libs)
      out="$out $libs"
    ;;
    --cxx-std)
      out="$out $cxx_std"
    ;;
    --prefix)
      out="$out $prefix"
    ;;
    --bindir)
      out="$out $bindir"
    ;;
    --libdir)
      out="$out $libdir"
    ;;
    --includedir)
      out="$out $includedir"
    ;;
    --use-root)
      out="$out $use_root"
    ;;
    --version)
      out="$out $version"
    ;;
    --git-revision)
      out="$out $git_revision"
    ;;

    # Flags that existed before the install layout was made relocatable.
    # Each gets its own message: "unknown argument" would be true but
    # useless, since the answer is different in every case.
    --incdir)
      echo "marley-config: --incdir was renamed --includedir" >&2
      exit 3
    ;;
    --topdir)
      echo "marley-config: --topdir was renamed --prefix. It reported" \
        "\$MARLEY, a source directory; it now reports the installation" \
        "prefix." >&2
      exit 3
    ;;
    --datadir)
      echo "marley-config: --datadir was removed. The data directory is" \
        "resolved at run time relative to the library, so baking it into a" \
        "build defeats that. Use 'marley info' to see the resolved value." >&2
      exit 3
    ;;
    --srcdir)
      echo "marley-config: --srcdir was removed. An installed MARLEY has no" \
        "source directory." >&2
      exit 3
    ;;

    *)
      echo "Unknown argument \"${arg}\"!" >&2
      echo "${usage_msg}" >&2
      exit 3
  esac
done

echo $out
