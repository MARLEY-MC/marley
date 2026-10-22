#!/bin/bash
#
# Throwaway test driver for the relocatable-install work. Not part of the
# build; delete before the clean replay.
#
#   ./check-install.sh
#
# Writes check-install.log (the interesting output) and
# check-install.build.log (compiler noise, consulted only on failure).

cd "$( dirname "${BASH_SOURCE[0]}" )" || exit 1
REPO="$PWD"
LOG="$REPO/check-install.log"
BUILDLOG="$REPO/check-install.build.log"

: > "$LOG"
: > "$BUILDLOG"

exec 3>&1                      # terminal

# Sections A6, A7 and B4 rename build/ out of the way to prove the install
# does not reach back into it. An interrupt there would leave the tree in
# that state, so restore it however this script exits.
restore_build() {
  if [ -d "$REPO/build.stashed" ] && [ ! -d "$REPO/build" ]; then
    mv "$REPO/build.stashed" "$REPO/build"
    printf 'restored build/ from build.stashed\n' >&3
  fi
}
trap restore_build EXIT INT TERM

log()  { printf '%s\n' "$*" >> "$LOG"; printf '%s\n' "$*" >&3; }
quiet(){ printf '%s\n' "$*" >> "$LOG"; }
head2(){ quiet ""; quiet "=== $* ==="; printf '>> %s\n' "$*" >&3; }

# Run a command, capturing both the command line and its output into the log.
run() {
  quiet ""
  quiet "\$ $*"
  "$@" >> "$LOG" 2>&1 < /dev/null
  local rc=$?
  [ $rc -ne 0 ] && quiet "[exit $rc]"
  return $rc
}

# A command that is SUPPOSED to fail. Without this, a deliberate refusal --
# "make uninstall" declining to proceed with no tty to confirm on -- is
# indistinguishable in the output from something actually breaking.
expect_fail() {
  quiet ""
  quiet "\$ $*   (expected to fail)"
  "$@" >> "$LOG" 2>&1 < /dev/null
  local rc=$?
  if [ $rc -eq 0 ]; then
    log "UNEXPECTED SUCCESS: $*"
  else
    quiet "[exit $rc, as expected]"
  fi
}

# Build quietly; on failure, splice the tail of the build log into the main
# log so the failure is visible without opening a second file.
build() {
  quiet ""
  quiet "\$ $*"
  echo "### $*" >> "$BUILDLOG"
  # The status is captured immediately. After a bare "if cmd; then ...; fi"
  # whose condition fails, $? is the *if statement's* status, which is 0 --
  # so reading it after the fi reported every failure as "[exit 0]".
  "$@" >> "$BUILDLOG" 2>&1 < /dev/null
  local rc=$?
  if [ $rc -eq 0 ]; then
    quiet "[ok]"
    return 0
  fi
  quiet "[exit $rc] -- last 40 lines of check-install.build.log:"
  tail -40 "$BUILDLOG" >> "$LOG"
  printf '!! FAILED: %s (see check-install.log)\n' "$*" >&3
  return $rc
}

rm -rf /tmp/m3 /tmp/m3-moved /tmp/m4 /tmp/m5 /tmp/marley-reloc
mkdir -p /tmp/marley-reloc

# A short job that still exercises the whole data-loading path: mass table,
# structure data, optical model, both .react families and their satellite
# files. Twenty events is enough to reach all of it.
cat > /tmp/marley-reloc/job.js <<'EOF'
{
  seed: 123456,
  target: { nuclides: [ 1000180400 ], atom_fractions: [ 1.0 ] },
  reactions: [
    "ve40ArCC_HF-CRPA.react",
    "ve40ArCC_Bhattacharya2009-Discrete.react",
    "ES.react"
  ],
  source: { type: "mono", neutrino: "ve", energy: 35.0 },
  direction: { x: 0.0, y: 0.0, z: 1.0 },
  generate: {
    events: 20,
    // "force" is the key that suppresses the overwrite prompt; "mode":
    // "overwrite" is the default and still prompts. Getting this wrong
    // hangs the run on a prompt nobody can see.
    output: [ { file: "events.hepmc3", format: "ascii", force: true } ],
  },
}
EOF

log "MARLEY install check -- $( date )"
log "repo: $REPO"
log "HEAD: $( git -C "$REPO" log --oneline -1 )"

# ---------------------------------------------------------------------------
# Configuration A: CMake, external HepMC3, external GSL, ROOT
# ---------------------------------------------------------------------------
head2 "A. CMake path -> /tmp/m3"

build make clean
build make -j10
build make install prefix=/tmp/m3

head2 "A1. installed tree (top two levels)"
run sh -c "find /tmp/m3 -maxdepth 2 -mindepth 1 -not -path '*/data/*' | sort"

head2 "A2. marley-config"
run /tmp/m3/bin/marley-config --cflags
run /tmp/m3/bin/marley-config --libs
run /tmp/m3/bin/marley-config --prefix
run /tmp/m3/bin/marley-config --bindir
run /tmp/m3/bin/marley-config --libdir
run /tmp/m3/bin/marley-config --includedir
run /tmp/m3/bin/marley-config --cxx-std
run /tmp/m3/bin/marley-config --version
run /tmp/m3/bin/marley-config --use-root
quiet ""
quiet "--- generated header, as installed ---"
sed -n '1,25p' /tmp/m3/bin/marley-config >> "$LOG" 2>&1
quiet ""
quiet "--- removed flags should each explain themselves ---"
expect_fail /tmp/m3/bin/marley-config --incdir
expect_fail /tmp/m3/bin/marley-config --datadir
expect_fail /tmp/m3/bin/marley-config --topdir
expect_fail /tmp/m3/bin/marley-config --srcdir

head2 "A3. marley.pc"
quiet ""
quiet "--- file ---"
cat /tmp/m3/lib/pkgconfig/marley.pc >> "$LOG" 2>&1
quiet ""
quiet "--- via pkg-config ---"
run env PKG_CONFIG_PATH=/tmp/m3/lib/pkgconfig pkg-config --cflags marley
run env PKG_CONFIG_PATH=/tmp/m3/lib/pkgconfig pkg-config --libs marley
run env PKG_CONFIG_PATH=/tmp/m3/lib/pkgconfig pkg-config --print-requires marley
run env PKG_CONFIG_PATH=/tmp/m3/lib/pkgconfig pkg-config --print-requires-private marley

head2 "A4. rootmap header keys (expect 7, all prefixed)"
run grep '^header' /tmp/m3/lib/libMARLEY.rootmap

head2 "A5. linkage of the installed artifacts"
if command -v otool > /dev/null; then
  run otool -L /tmp/m3/bin/marley
  run otool -l /tmp/m3/bin/marley
  run otool -L /tmp/m3/lib/libMARLEY.dylib
else
  run ldd /tmp/m3/bin/marley
  run readelf -d /tmp/m3/bin/marley
  run readelf -d /tmp/m3/lib/libMARLEY.so
fi

head2 "A6. relocation: run with the build tree moved away"
mv "$REPO/build" "$REPO/build.stashed"
mkdir -p /tmp/marley-reloc/a6
( cd /tmp/marley-reloc/a6 \
  && /tmp/m3/bin/marley ../job.js ) >> "$LOG" 2>&1 < /dev/null
quiet "[exit $?]"
mv "$REPO/build.stashed" "$REPO/build"

head2 "A7. same run, from a relocated copy of the whole prefix"
cp -R /tmp/m3 /tmp/m3-moved
mv "$REPO/build" "$REPO/build.stashed"
mkdir -p /tmp/marley-reloc/a7
( cd /tmp/marley-reloc/a7 \
  && /tmp/m3-moved/bin/marley ../job.js ) >> "$LOG" 2>&1 < /dev/null
quiet "[exit $?]"
run /tmp/m3-moved/bin/marley-config --libdir
mv "$REPO/build.stashed" "$REPO/build"

# ---------------------------------------------------------------------------
# Configuration B: GNU Make, vendored HepMC3
# ---------------------------------------------------------------------------
head2 "B. Make path, vendored HepMC3 -> /tmp/m4"

build make clean
build make IGNORE_CMAKE=1 IGNORE_HEPMC3=1 -j10
build make IGNORE_CMAKE=1 IGNORE_HEPMC3=1 install prefix=/tmp/m4

head2 "B1. installed tree (top two levels)"
run sh -c "find /tmp/m4 -maxdepth 2 -mindepth 1 -not -path '*/data/*' | sort"
run ls /tmp/m4/include/marley/builtin

head2 "B2. marley-config (vendored branch of every conditional)"
run /tmp/m4/bin/marley-config --cflags
run /tmp/m4/bin/marley-config --libs
run /tmp/m4/bin/marley-config --prefix
run /tmp/m4/bin/marley-config --cxx-std

head2 "B3. marley.pc (vendored: no Requires, -lMARLEYHepMC3)"
cat /tmp/m4/lib/pkgconfig/marley.pc >> "$LOG" 2>&1
run env PKG_CONFIG_PATH=/tmp/m4/lib/pkgconfig pkg-config --cflags marley
run env PKG_CONFIG_PATH=/tmp/m4/lib/pkgconfig pkg-config --libs marley

head2 "B3a. build-tree copy of the scripts, vendored build"
quiet "builtin_rel must be include/builtin here (install tree uses"
quiet "include/marley/builtin), and hepmc3_include must be empty."
run sh -c "grep -E \"(builtin_rel|hepmc3_include)=\" \
  build/bin/setup_marley.sh build/bin/mroot"

head2 "B4. relocation"
mv "$REPO/build" "$REPO/build.stashed"
mkdir -p /tmp/marley-reloc/b4
( cd /tmp/marley-reloc/b4 \
  && /tmp/m4/bin/marley ../job.js ) >> "$LOG" 2>&1 < /dev/null
quiet "[exit $?]"
mv "$REPO/build.stashed" "$REPO/build"

head2 "D. Make path, EXTERNAL HepMC3 -> /tmp/m5"
# The combination neither A nor B covers. It is the only configuration in
# which -DMARLEY_FOUND_HEPMC3 can reach the published flags: the CMake
# config-package path never appends it, and a vendored build has no external
# HepMC3 flags at all. Neither must appear in --cflags (design spec 6.2).

build make clean
build make IGNORE_CMAKE=1 -j10
build make IGNORE_CMAKE=1 install prefix=/tmp/m5

run /tmp/m5/bin/marley-config --cflags
run /tmp/m5/bin/marley-config --libs
quiet ""
quiet "--- these two must NOT appear above ---"
if /tmp/m5/bin/marley-config --cflags --libs \
     | grep -q -e -DMARLEY_FOUND_HEPMC3 -e -DMARLEY_FOUND_GSL; then
  log "LEAK: a MARLEY_FOUND_* build switch reached the published flags"
else
  quiet "clean: neither MARLEY_FOUND_HEPMC3 nor MARLEY_FOUND_GSL is published"
fi
cat /tmp/m5/lib/pkgconfig/marley.pc >> "$LOG" 2>&1
run env PKG_CONFIG_PATH=/tmp/m5/lib/pkgconfig pkg-config --cflags marley
run env PKG_CONFIG_PATH=/tmp/m5/lib/pkgconfig pkg-config --libs marley

head2 "G. generated shell scripts (design spec 6.4)"
quiet ""
quiet "--- installed mroot and setup_marley.sh, external-HepMC3 build ---"
cat /tmp/m5/bin/setup_marley.sh >> "$LOG" 2>&1
quiet ""
run sh -c "grep -E \"^(exec|lib_rel|include_rel|builtin_rel)\" /tmp/m5/bin/mroot"
quiet ""
quiet "--- vendored build: builtin_rel must be include/marley/builtin ---"
run sh -c "grep -E \"builtin_rel=\" /tmp/m4/bin/mroot /tmp/m4/bin/setup_marley.sh"
quiet ""
quiet "--- external-HepMC3 build: hepmc3_include must be the HepMC3 prefix ---"
run sh -c "grep -E \"(builtin_rel|hepmc3_include)=\" \
  /tmp/m5/bin/setup_marley.sh /tmp/m5/bin/mroot"
quiet ""
quiet "--- sourcing the installed script must set no MARLEY variable ---"
run sh -c "unset MARLEY; . /tmp/m5/bin/setup_marley.sh >/dev/null 2>&1; \
  echo \"MARLEY=[\${MARLEY-<unset>}]\"; \
  echo \"PATH head=\$( echo \"\$PATH\" | cut -d: -f1 )\"; \
  echo \"ROOT_INCLUDE_PATH head=\$( echo \"\$ROOT_INCLUDE_PATH\" | cut -d: -f1 )\""
quiet ""
quiet "--- repo-root wrapper, with no build directory present ---"
mv "$REPO/build" "$REPO/build.stashed"
run sh -c ". $REPO/setup_marley.sh"
mv "$REPO/build.stashed" "$REPO/build"
quiet ""
quiet "--- mroot must load the library and parse a MARLEY header ---"
for p in /tmp/m5 /tmp/m4; do
  quiet ""
  quiet "\$ $p/bin/mroot -b -q  (parse marley/EventFileReader.hh)"
  # The USE_ROOT probe must come first: Declare() returns false if the
  # #error fires, so a 0 means cling disagrees with the compiled library
  # about EventFileReader's layout.
  "$p/bin/mroot" -b -q \
    -e 'printf("USE_ROOT reached cling: %d\n", (int)gInterpreter->Declare("#ifndef USE_ROOT\n#error USE_ROOT is not defined\n#endif\n"));' \
    -e '#include "marley/EventFileReader.hh"' \
    -e 'printf("HEADER PARSED OK\n");' \
    -e '(int)0;' >> "$LOG" 2>&1 < /dev/null
  quiet "[exit $?]"
  # The trailing '(int)0' matters: root -q exits with the value of the LAST
  # -e expression, and printf returns its character count -- which made a
  # clean run look like it exited 17.
done

head2 "H. marley info (design spec 2.5, 6.1)"
quiet ""
quiet "--- plain: library, roots with provenance, cwd last ---"
run /tmp/m5/bin/marley info
quiet ""
quiet "--- a name that resolves, with the subdir its real caller passes ---"
run /tmp/m5/bin/marley info --resolve --subdir react ES.react
quiet ""
quiet "--- the same name WITHOUT --subdir: must fail and suggest ---"
expect_fail /tmp/m5/bin/marley info --resolve ES.react
quiet ""
quiet "--- a name that sits at the top of the data root ---"
run /tmp/m5/bin/marley info --resolve mass_table.js
quiet ""
quiet "--- a name that exists nowhere ---"
expect_fail /tmp/m5/bin/marley info --resolve no_such_file.js
quiet ""
quiet "--- MARLEY_SEARCH_PATH roots are flagged, and a bad one is marked ---"
run env MARLEY_SEARCH_PATH=/tmp/marley-reloc:/tmp/definitely-not-here \
  /tmp/m5/bin/marley info
quiet ""
quiet "--- from the relocated prefix: roots must follow the copy ---"
run /tmp/m3-moved/bin/marley info
quiet ""
quiet "--- argument errors ---"
expect_fail /tmp/m5/bin/marley info --resolve
expect_fail /tmp/m5/bin/marley info ES.react
expect_fail /tmp/m5/bin/marley info --subdir react
expect_fail /tmp/m5/bin/marley info --bogus
quiet ""
quiet "--- help, and presence in the top-level command list ---"
run /tmp/m5/bin/marley info --help
run sh -c "/tmp/m5/bin/marley help | grep -A1 info"

head2 "I. the C++ standard must be one value (design spec 6.2)"
quiet ""
quiet "Both build systems saw the same ROOT, so --cxx-std must agree. It did"
quiet "not before: CMake reported c++17 while Make reported c++20, because"
quiet "ROOT's PUBLIC compile features raised the standard behind the value"
quiet "CMake was reporting."
for p in /tmp/m3 /tmp/m4 /tmp/m5; do
  quiet ""
  quiet "\$ $p/bin/marley-config --cxx-std"
  "$p/bin/marley-config" --cxx-std >> "$LOG" 2>&1 < /dev/null
done
quiet ""
quiet "--- and marley.pc must say the same thing ---"
run sh -c "grep -h '^Cflags:' /tmp/m3/lib/pkgconfig/marley.pc \
  /tmp/m4/lib/pkgconfig/marley.pc /tmp/m5/lib/pkgconfig/marley.pc \
  | sed 's/.*\\(-std=[^ ]*\\).*/\\1/'"
quiet ""
quiet "--- exactly one -std= in the published cflags, not two ---"
run sh -c "/tmp/m5/bin/marley-config --cflags | tr ' ' '\\n' | grep -c -- -std="
quiet ""
quiet "--- an explicitly requested standard must not be silently downgraded ---"
quiet "(c++23 against a C++20 ROOT: the old code fired the c++20 block and"
quiet "overrode the request back down)"
build make clean
build make IGNORE_CMAKE=1 CXX_STD=c++23 -j10
run sh -c "./build/bin/marley-config --cxx-std"

head2 "J. the CMake package must work for a consumer"
quiet ""
quiet "A find_package(MARLEY) that configures, builds and links is the only"
quiet "real check: the config file previously had its find_dependency calls"
quiet "commented out, so it resolved and failed later at link time."
quiet ""
quiet "--- what the installed package actually contains ---"
run sh -c "ls /tmp/m3/lib/cmake/MARLEY"
run sh -c "cat /tmp/m3/lib/cmake/MARLEY/MARLEYConfig.cmake"
run sh -c "grep -nE 'add_library|_IMPORT_CHECK_TARGETS|set_target_properties' \
  /tmp/m3/lib/cmake/MARLEY/MARLEYTargets.cmake | head -20"

rm -rf /tmp/marley-consumer && mkdir -p /tmp/marley-consumer
cat > /tmp/marley-consumer/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.18)
project(marley_consumer LANGUAGES CXX)
find_package(MARLEY REQUIRED)
# Report every imported target the package defined, so that a wrong guess
# at the name is distinguishable from a package that exported nothing.
get_property(_tgts DIRECTORY PROPERTY IMPORTED_TARGETS)
message(STATUS "imported targets after find_package(MARLEY): ${_tgts}")
message(STATUS "MARLEY_VERSION=${MARLEY_VERSION} MARLEY_FOUND=${MARLEY_FOUND}")
add_executable(consumer main.cc)
target_link_libraries(consumer PRIVATE MARLEY::marley_lib)
EOF
cat > /tmp/marley-consumer/main.cc <<'EOF'
#include "marley/FileManager.hh"
#include <iostream>
int main() {
  const auto& fm = marley::FileManager::Instance();
  std::cout << "consumer sees " << fm.roots().size() << " data root(s)\n";
  return 0;
}
EOF
quiet ""
quiet "--- configure, build, run ---"
run cmake -S /tmp/marley-consumer -B /tmp/marley-consumer/build \
  -DCMAKE_PREFIX_PATH=/tmp/m3
run cmake --build /tmp/marley-consumer/build
run /tmp/marley-consumer/build/consumer
quiet ""
quiet "--- the consumer must inherit the C++ standard from the target ---"
run sh -c "grep -o -- '-std=[^ ]*' \
  /tmp/marley-consumer/build/CMakeFiles/consumer.dir/flags.make | head -2"

head2 "C. no build-tree or source-tree path in any published interface"
quiet ""
quiet "Searching the installed marley-config and marley.pc for \"$REPO\"."
quiet "Any hit is a bug: an installed MARLEY must not reference the tree it"
quiet "was built from."
for f in /tmp/m3/bin/marley-config /tmp/m3/lib/pkgconfig/marley.pc \
         /tmp/m4/bin/marley-config /tmp/m4/lib/pkgconfig/marley.pc \
         /tmp/m5/bin/marley-config /tmp/m5/lib/pkgconfig/marley.pc; do
  if grep -n -- "$REPO" "$f" >> "$LOG" 2>&1; then
    log "LEAK: $f references the build tree"
  else
    quiet "clean: $f"
  fi
done
quiet ""
quiet "--- and the emitted flags, which is what consumers actually see ---"
for p in /tmp/m3 /tmp/m4 /tmp/m5; do
  out="$( "$p/bin/marley-config" --cflags --libs 2>&1 )"
  quiet "$p: $out"
  case "$out" in
    *"$REPO"*) log "LEAK: $p/bin/marley-config emits a build-tree path" ;;
    *)         quiet "clean: $p" ;;
  esac
done

head2 "F. no rpath anywhere in any marley.pc"
quiet ""
quiet "marley-config --libs may carry an rpath; marley.pc must not, and that"
quiet "includes the one root-config emits for ROOT's own library directory."
# Only the Cflags:/Libs: directives, not the whole file -- the template's
# own comment explains the rpath rule and contains the string.
for f in /tmp/m3/lib/pkgconfig/marley.pc /tmp/m4/lib/pkgconfig/marley.pc \
         /tmp/m5/lib/pkgconfig/marley.pc; do
  if grep -E '^(Cflags|Libs):' "$f" | grep -n -- "-Wl,-rpath" >> "$LOG" 2>&1
  then
    log "RPATH LEAK: $f"
  else
    quiet "clean: $f"
  fi
done

head2 "E. uninstall by receipt (design spec 7)"
quiet ""
quiet "--- the receipts both installers wrote ---"
run sh -c "cat /tmp/m3/share/marley/install-receipt"
run sh -c "cat /tmp/m4/share/marley/install-receipt"
quiet ""
quiet "--- a user file placed in the prefix must survive ---"
mkdir -p /tmp/m4/bin/mine && : > /tmp/m4/bin/mine/keep.txt
quiet ""
quiet "--- make uninstall prompts, and a non-tty answer aborts safely ---"
expect_fail make IGNORE_CMAKE=1 IGNORE_HEPMC3=1 uninstall prefix=/tmp/m4
run find /tmp/m4 -maxdepth 2
quiet ""
quiet "--- the same prefix, removed non-interactively as CI would ---"
run src/scripts/uninstall-marley --prefix /tmp/m4 --yes
run find /tmp/m4
quiet ""
quiet "--- CMake-path install, removed by the same script ---"
run src/scripts/uninstall-marley --prefix /tmp/m3 --yes
run find /tmp/m3
quiet ""
quiet "--- and the refusals ---"
expect_fail src/scripts/uninstall-marley
expect_fail src/scripts/uninstall-marley --prefix /
expect_fail src/scripts/uninstall-marley --prefix /tmp/m3 --yes

log ""
log "Done. Wrote check-install.log (and check-install.build.log)."
