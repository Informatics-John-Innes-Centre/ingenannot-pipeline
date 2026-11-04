#! /usr/bin/env bash
# strip_utrs.sh - remove UTRs from a GFF3 and resize features to match.
#
# For every transcript with CDS children:
#   * UTR features (five_prime_UTR, three_prime_UTR, UTR, ...) are dropped
#   * transcript start/end -> min/max of its CDS
#   * exons are clipped to that CDS span (exons entirely outside are dropped)
# Genes are resized to the union of their remaining transcripts.
# Non-coding transcripts are kept unchanged unless --drop-noncoding is given.
# Exons with several parents (Parent=a,b) are split into one line per parent.
#
# Usage:
#   strip_utrs.sh in.gff3 > out.gff3
#   strip_utrs.sh --drop-noncoding in.gff3 > out.gff3

set -euo pipefail

drop=0
gff=""
for arg in "$@"; do
    case "$arg" in
        --drop-noncoding) drop=1 ;;
        -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
        *) gff="$arg" ;;
    esac
done

if [[ -z "$gff" || ! -f "$gff" ]]; then
    echo "Usage: $0 [--drop-noncoding] in.gff3 > out.gff3" >&2
    exit 1
fi

# The file is read twice: pass 1 collects spans, pass 2 writes output.
awk -F'\t' -v OFS='\t' -v drop="$drop" '
function attr(s, key,   n, i, kv, parts) {
    n = split(s, parts, ";")
    for (i = 1; i <= n; i++) {
        if (index(parts[i], key "=") == 1) return substr(parts[i], length(key) + 2)
    }
    return ""
}
# rebuild attribute string with a single Parent (and suffixed ID if split)
function reattr(s, p, multi,   n, i, parts, out, kv) {
    n = split(s, parts, ";"); out = ""
    for (i = 1; i <= n; i++) {
        if (parts[i] == "") continue
        kv = parts[i]
        if (index(kv, "Parent=") == 1) kv = "Parent=" p
        else if (multi && index(kv, "ID=") == 1) kv = kv "." p
        out = out (out == "" ? "" : ";") kv
    }
    return out
}
function is_utr(t) { return tolower(t) ~ /utr/ }

# finish pass 1: transcript spans -> gene spans
function summarise(   id, n, i, ps, s, e) {
    for (id in is_tx) {
        if (id in cmin) { s = cmin[id]; e = cmax[id] }
        else if (drop) { dropped[id] = 1; n_nc++; continue }
        else { s = st[id]; e = en[id] }
        n = split(par[id], ps, ",")
        for (i = 1; i <= n; i++) {
            if (!(ps[i] in gmin) || s < gmin[ps[i]]) gmin[ps[i]] = s
            if (!(ps[i] in gmax) || e > gmax[ps[i]]) gmax[ps[i]] = e
        }
    }
    for (id in is_tx) {
        n = split(par[id], ps, ",")
        for (i = 1; i <= n; i++) if (!(ps[i] in gmin)) empty[ps[i]] = 1
    }
    done = 1
}

# ---------- pass 1 ----------
FNR == NR {
    if ($0 ~ /^#/ || NF != 9) next
    id = attr($9, "ID"); p = attr($9, "Parent")
    if (id != "") { st[id] = $4; en[id] = $5; par[id] = p }
    if (p == "") next
    n = split(p, ps, ",")
    for (i = 1; i <= n; i++) {
        if ($3 == "exon" || $3 == "CDS") is_tx[ps[i]] = 1
        if ($3 == "CDS") {
            if (!(ps[i] in cmin) || $4 + 0 < cmin[ps[i]]) cmin[ps[i]] = $4 + 0
            if (!(ps[i] in cmax) || $5 + 0 > cmax[ps[i]]) cmax[ps[i]] = $5 + 0
        }
    }
    next
}

# ---------- pass 2 ----------
!done { summarise() }

$0 ~ /^#/ || NF != 9 { print; next }

{
    id = attr($9, "ID"); p = attr($9, "Parent")
    if (is_utr($3)) { n_utr++; next }
    if ((id != "" && (id in dropped || id in empty))) next

    np = (p == "") ? 0 : split(p, ps, ",")
    if (np > 0) {
        alldrop = 1
        for (i = 1; i <= np; i++) if (!(ps[i] in dropped)) alldrop = 0
        if (alldrop) next
    }

    if ($3 == "exon" && np > 0) {
        s0 = $4 + 0; e0 = $5 + 0; a0 = $9
        for (i = 1; i <= np; i++) {
            q = ps[i]
            if (q in dropped) continue
            s = s0; e = e0
            if (q in cmin) {
                if (cmin[q] > s) s = cmin[q]
                if (cmax[q] < e) e = cmax[q]
                if (s > e) { n_ex++; continue }
            }
            $4 = s; $5 = e
            $9 = (np > 1) ? reattr(a0, q, 1) : a0
            print
        }
        next
    }

    if (id in is_tx) {
        if (id in cmin) { $4 = cmin[id]; $5 = cmax[id] }
    } else if (id in gmin) {
        $4 = gmin[id]; $5 = gmax[id]
    }
    print
}

END {
    printf "dropped %d UTR features, %d exons outside CDS, %d non-coding transcripts\n",
        n_utr, n_ex, n_nc > "/dev/stderr"
}
' "$gff" "$gff"
