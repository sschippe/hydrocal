#!/bin/bash
# SPDX-License-Identifier: MIT
# Run a single test case. Args: <hcin-file> <result-file>
# Writes "ok" or "FAILED (reason)" to result-file.

hcin=$(readlink -f "$1")
result_file="$2"
name=$(basename "$hcin" .hcin)
refdir=$(dirname "$hcin")
# Use an absolute logdir so the output checks below still resolve correctly
# after we cd into $srcdir (previously logdir was relative and the checks
# looked in the wrong directory, making every test report "missing").
logdir="$(pwd)/test"
srcdir=$(dirname "$(readlink -f "$0")")

mkdir -p "$logdir/$name"

# copy support files from the examples directory
for f in Ne2deconv.inp O2DN1.bdr Xe53DRtheo4.dat; do
  [ -f "$refdir/$f" ] && cp "$refdir/$f" "$logdir/$name/"
done

# find and copy hydrocal binary
bp=""
for b in "$srcdir/hydrocal.exe" "$srcdir/hydrocal" "$srcdir/build/hydrocal"; do
  [ -f "$b" ] && bp="$b" && break
done
if [ -z "$bp" ]; then
  echo "FAILED (hydrocal binary not found)" > "$result_file"
  exit 1
fi
cp -f "$bp" "$logdir/$name/hydrocal"

# Ds3p32DiracRR validates the DIR dipole cross section against its
# reference; the fully retarded default mode is too slow for Z=110 and
# would change the expected values.  Force the dipole path for this case.
case "$name" in
  Ds3p32DiracRR) export HYDROCAL_DIRAC_RETARDED=0 ;;
  *) unset HYDROCAL_DIRAC_RETARDED HYDROCAL_DIRAC_RETARDED_NUMERIC ;;
esac

# run hydrocal from the test directory
cd "$logdir/$name"
./hydrocal < "$hcin" > hydrocal.log 2>&1
rc=$?
cd "$srcdir"

if [ $rc -ne 0 -a $rc -ne 1 ]; then
  echo "FAILED (exit code $rc)" > "$result_file"
  exit 1
fi

case "$name" in
  S3soft)
    for out in S3soft.fn S3soft.fnl; do
      if [ ! -f "$logdir/$name/$out" ]; then
        echo "FAILED (missing $out)" > "$result_file"; exit 1
      fi
      grep -v '^#' "$logdir/$name/$out" > /tmp/s3soft_new_$$
      grep -v '^#' "$refdir/$out" > /tmp/s3soft_ref_$$
      diff -q /tmp/s3soft_new_$$ /tmp/s3soft_ref_$$ >/dev/null 2>&1
      rc=$?
      rm -f /tmp/s3soft_new_$$ /tmp/s3soft_ref_$$
      if [ $rc -ne 0 ]; then
        echo "FAILED ($out differs)" > "$result_file"; exit 1
      fi
    done
    echo "ok" > "$result_file"
    ;;
  U92RR)
    if [ ! -f "$logdir/$name/U92RR.crr" ]; then
      echo "FAILED (missing U92RR.crr)" > "$result_file"; exit 1
    fi
    grep -v '^#' "$logdir/$name/U92RR.crr" > "$logdir/$name/U92RR.crr.clean"
    grep -v '^#' "$refdir/U92RR.crr" > "$logdir/$name/U92RR.crr.ref"
    diff -q "$logdir/$name/U92RR.crr.clean" "$logdir/$name/U92RR.crr.ref" >/dev/null 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  CG|ThreeJ|SixJ|NineJ|OsciBB|OsciBC)
    case "$name" in
      CG)      outpat="Clebsch Gordan.*=";    ctx=("-B0" "-A0");;
      ThreeJ)  outpat="3-J Symbol.*=";        ctx=("-B1" "-A1");;
      SixJ)    outpat="6-J Symbol.*=";        ctx=("-B1" "-A1");;
      NineJ)   outpat="9-J Symbol.*=";        ctx=("-B1" "-A1");;
      OsciBB)  outpat="fBB =";                ctx=("-B0" "-A0");;
      OsciBC)  outpat="nfBC =";               ctx=("-B0" "-A0");;
    esac
    grep "${ctx[@]}" "$outpat" "$logdir/$name/hydrocal.log" | sed 's/[[:space:]]*$//' > /tmp/${name}_new_$$
    sed 's/[[:space:]]*$//' "$refdir/$name.ref" > /tmp/${name}_ref_$$
    if [ ! -s /tmp/${name}_new_$$ ]; then
      echo "FAILED (no output for '$pat')" > "$result_file"
      rm -f /tmp/${name}_new_$$ /tmp/${name}_ref_$$
      exit 1
    fi
    diff -q /tmp/${name}_new_$$ /tmp/${name}_ref_$$ >/dev/null 2>&1
    rc=$?
    rm -f /tmp/${name}_new_$$ /tmp/${name}_ref_$$
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (result differs)" > "$result_file"
    fi
    ;;
  Hlifetimes)
    if [ ! -f "$logdir/$name/Hlifetimes.tau" ]; then
      echo "FAILED (missing Hlifetimes.tau)" > "$result_file"; exit 1
    fi
    diff -q "$logdir/$name/Hlifetimes.tau" "$refdir/Hlifetimes.tau" >/dev/null 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  U91DiracRates)
    if [ ! -f "$logdir/$name/U91DiracRates.dtr" ]; then
      echo "FAILED (missing U91DiracRates.dtr)" > "$result_file"; exit 1
    fi
    grep -v '^#' "$logdir/$name/U91DiracRates.dtr" > "$logdir/$name/U91DiracRates.dtr.clean"
    grep -v '^#' "$refdir/U91DiracRates.dtr" > "$logdir/$name/U91DiracRates.dtr.ref"
    diff -q "$logdir/$name/U91DiracRates.dtr.clean" "$logdir/$name/U91DiracRates.dtr.ref" >/dev/null 2>&1
    rc=$?
    rm -f "$logdir/$name/U91DiracRates.dtr.clean" "$logdir/$name/U91DiracRates.dtr.ref"
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  Ds3p32DiracRR)
    if [ ! -f "$logdir/$name/Ds3p32DiracRR.sig" ]; then
      echo "FAILED (missing Ds3p32DiracRR.sig)" > "$result_file"; exit 1
    fi
    diff -q "$logdir/$name/Ds3p32DiracRR.sig" "$refdir/Ds3p32DiracRR.sig" >/dev/null 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  O2cooler)
    if [ ! -f "$logdir/$name/O2DN11.mcr" ]; then
      echo "FAILED (missing O2DN11.mcr)" > "$result_file"; exit 1
    fi
    lines=$(wc -l < "$logdir/$name/O2DN11.mcr")
    if [ $lines -lt 50 ]; then
      echo "FAILED (too short: $lines lines)" > "$result_file"
    else
      echo "ok ($lines lines)" > "$result_file"
    fi
    ;;
  Xe53lorentz)
    if [ ! -f "$logdir/$name/Xe53lorentz.cdr" ]; then
      echo "FAILED (missing Xe53lorentz.cdr)" > "$result_file"; exit 1
    fi
    grep -v '^#' "$logdir/$name/Xe53lorentz.cdr" > "$logdir/$name/Xe53lorentz.cdr.clean"
    grep -v '^#' "$refdir/Xe53lorentz.cdr" > "$logdir/$name/Xe53lorentz.cdr.ref"
    diff -q "$logdir/$name/Xe53lorentz.cdr.clean" "$logdir/$name/Xe53lorentz.cdr.ref" >/dev/null 2>&1
    rc=$?
    rm -f "$logdir/$name/Xe53lorentz.cdr.clean" "$logdir/$name/Xe53lorentz.cdr.ref"
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  Pb82RRphotons)
    for out in Pb82RRphotons.RRlines Pb82RRphotons.dat; do
      if [ ! -f "$logdir/$name/$out" ]; then
        echo "FAILED (missing $out)" > "$result_file"; exit 1
      fi
      grep -v '^#' "$logdir/$name/$out" > /tmp/pb82rr_new_$$
      grep -v '^#' "$refdir/$out" > /tmp/pb82rr_ref_$$
      diff -q /tmp/pb82rr_new_$$ /tmp/pb82rr_ref_$$ >/dev/null 2>&1
      rc=$?
      rm -f /tmp/pb82rr_new_$$ /tmp/pb82rr_ref_$$
      if [ $rc -ne 0 ]; then
        echo "FAILED ($out differs)" > "$result_file"; exit 1
      fi
    done
    echo "ok" > "$result_file"
    ;;
  Ne2deconv)
    if [ ! -f "$logdir/$name/Ne2deconv.dcv" ]; then
      echo "FAILED (missing Ne2deconv.dcv)" > "$result_file"; exit 1
    fi
    grep -v '^#' "$logdir/$name/Ne2deconv.dcv" > "$logdir/$name/Ne2deconv.dcv.clean"
    grep -v '^#' "$refdir/Ne2deconv.dcv" > "$logdir/$name/Ne2deconv.dcv.ref"
    diff -q "$logdir/$name/Ne2deconv.dcv.clean" "$logdir/$name/Ne2deconv.dcv.ref" >/dev/null 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  CNplus)
    if [ ! -f "$logdir/$name/CNplus.fcf" ]; then
      echo "FAILED (missing CNplus.fcf)" > "$result_file"; exit 1
    fi
    grep -v '^#' "$logdir/$name/CNplus.fcf" > "$logdir/$name/CNplus.fcf.clean"
    grep -v '^#' "$refdir/CNplus.fcf" > "$logdir/$name/CNplus.fcf.ref"
    diff -q "$logdir/$name/CNplus.fcf.clean" "$logdir/$name/CNplus.fcf.ref" >/dev/null 2>&1
    rc=$?
    rm -f "$logdir/$name/CNplus.fcf.clean" "$logdir/$name/CNplus.fcf.ref"
    if [ $rc -eq 0 ]; then echo "ok" > "$result_file"
    else echo "FAILED (data differs)" > "$result_file"
    fi
    ;;
  *)
    echo "FAILED (unknown test)" > "$result_file"
    ;;
esac

# The success branches above fall through to here, so translate the verdict
# written to the result file into the exit status.  Without this every test
# reported success to CTest regardless of whether the numbers matched.
if grep -q '^FAILED' "$result_file"; then
  exit 1
fi
exit 0
