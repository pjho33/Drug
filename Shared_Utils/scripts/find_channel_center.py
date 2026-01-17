#!/usr/bin/env python3
import math
import sys


def _parse_pdb_coords(pdb_path, record_prefix, atom_name=None, resname=None):
    out = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith(record_prefix):
                continue
            if resname is not None and line[17:20].strip() != resname:
                continue
            if atom_name is not None and line[12:16].strip() != atom_name:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain_id = (line[21].strip() or '_')
            out.append((chain_id, x, y, z))
    return out


def _centroid(points):
    if not points:
        return None
    sx = sy = sz = 0.0
    for x, y, z in points:
        sx += x
        sy += y
        sz += z
    n = float(len(points))
    return (sx / n, sy / n, sz / n)


def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def main():
    pdb_path = sys.argv[1] if len(sys.argv) >= 2 else '/home/pjho3/projects/Drug/raw_data/4PYP_trimer.pdb'

    bng_atoms = _parse_pdb_coords(pdb_path, 'HETATM', resname='BNG')
    bng_by_chain = {}
    for chain_id, x, y, z in bng_atoms:
        bng_by_chain.setdefault(chain_id, []).append((x, y, z))

    if not bng_by_chain:
        print('No BNG found. Provide a PDB that includes bound glucose-analog (e.g., BNG) to use as pocket proxy.')
        return 2

    centroids = {ch: _centroid(coords) for ch, coords in bng_by_chain.items()}
    chains = sorted(centroids.keys())

    print(f'PDB: {pdb_path}')
    print('BNG pocket-proxy centroids:')
    for ch in chains:
        c = centroids[ch]
        print(f'  Chain {ch}: ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f})')

    print('Pairwise distances (BNG centroid):')
    for i in range(len(chains)):
        for j in range(i + 1, len(chains)):
            a = chains[i]
            b = chains[j]
            print(f'  {a}-{b}: {_dist(centroids[a], centroids[b]):.2f} Å')

    ca_atoms = _parse_pdb_coords(pdb_path, 'ATOM  ', atom_name='CA')
    ca_by_chain = {}
    for chain_id, x, y, z in ca_atoms:
        ca_by_chain.setdefault(chain_id, []).append((x, y, z))
    ca_centroids = {ch: _centroid(coords) for ch, coords in ca_by_chain.items()}

    print('Pocket direction proxy (BNG centroid - CA centroid):')
    for ch in chains:
        if ch not in ca_centroids or ca_centroids[ch] is None:
            continue
        b = centroids[ch]
        c = ca_centroids[ch]
        v = (b[0] - c[0], b[1] - c[1], b[2] - c[2])
        print(f'  Chain {ch}: ({v[0]:.2f}, {v[1]:.2f}, {v[2]:.2f}) |v|={_dist((0.0,0.0,0.0), v):.2f} Å')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
