#!/usr/bin/env python3
"""Copy this ALPGEN six-parton LHE and record its actual ME alpha_s in AQCDUP."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re


def number(text):
    return float(text.replace('D', 'E').replace('d', 'e'))


def convert(source, destination, lhapdf_base):
    if source.resolve() == destination.resolve():
        raise ValueError('The original LHE must not be overwritten')
    # These checks restrict the scale interpretation to the verified ALPGEN setup.
    header = []
    with source.open() as stream:
        for line in stream:
            header.append(line)
            if '<init>' in line:
                break
        else:
            raise ValueError('Missing LHE init block')
    parameters = {}
    for line in header:
        match = re.match(r'^\s*\d+\s+(\S+)\s+!\s*(\w+)\s*$', line)
        if match:
            parameters[match[2]] = number(match[1])
    expected = dict(ndns=315200, iqopt=1, qfac=1, njets=6, ickkw=0, ebeam=6800)
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise ValueError(f'Unsupported ALPGEN {key}: {parameters.get(key)}, expected {value}')

    import ROOT
    ROOT.gInterpreter.AddIncludePath(str(lhapdf_base / 'include'))
    if ROOT.gSystem.Load('libLHAPDF') < 0:
        raise RuntimeError('Cannot load CMSSW libLHAPDF')
    if not ROOT.gInterpreter.Declare('#include <LHAPDF/LHAPDF.h>'):
        raise RuntimeError('Cannot load LHAPDF headers')
    pdf = ROOT.LHAPDF.mkPDF('NNPDF31_lo_as_0130', 0)
    count = replaced = 0
    as_min, as_max = math.inf, -math.inf
    original_hash, corrected_hash = hashlib.sha256(), hashlib.sha256()
    temp = destination.with_name(destination.name + '.tmp')
    closed = False
    try:
        with source.open() as src, temp.open('w') as dst:
            def write(line):
                dst.write(line)
                corrected_hash.update(line.encode())

            for line in src:
                original_hash.update(line.encode())
                write(line)
                if '</LesHouchesEvents>' in line:
                    closed = True
                if line.strip() != '<event>':
                    continue
                block = []
                for line in src:
                    original_hash.update(line.encode())
                    block.append(line)
                    if line.strip() == '</event>':
                        break
                else:
                    raise ValueError('Unclosed LHE event')
                idx = next(i for i, line in enumerate(block) if line.strip())
                fields = block[idx].split()
                if len(fields) != 6:
                    raise ValueError('Invalid LHE event header')
                nup, scale, previous = int(fields[0]), number(fields[3]), number(fields[5])
                if not math.isfinite(scale) or scale <= 0:
                    raise ValueError('Invalid SCALUP')
                particles = [line.split() for line in block[idx + 1:idx + 1 + nup]]
                if len(particles) != nup or any(len(p) != 13 for p in particles):
                    raise ValueError('Invalid LHE particles')
                final = [p for p in particles if int(p[1]) == 1]
                if nup != 8 or len(final) != 6 or any(abs(int(p[0])) not in (1,2,3,4,5,21) for p in final):
                    raise ValueError('Expected exactly six QCD final-state partons')
                expected_scale = math.sqrt(sum(number(p[6])**2 + number(p[7])**2 for p in final) / 6)
                if not math.isclose(scale, expected_scale, rel_tol=5e-5):
                    raise ValueError(f'SCALUP {scale} disagrees with ME scale {expected_scale}')
                alpha = float(pdf.alphasQ(scale))
                if not math.isfinite(alpha) or not 0 < alpha < 1:
                    raise ValueError('Invalid LHAPDF alpha_s')
                if previous < 0:
                    token = list(re.finditer(r'\S+', block[idx]))[-1]
                    block[idx] = block[idx][:token.start()] + f'{alpha:.12e}' + block[idx][token.end():]
                    replaced += 1
                elif not math.isclose(previous, alpha, rel_tol=5e-5):
                    raise ValueError(f'Existing AQCDUP {previous} disagrees with ME alpha_s {alpha}')
                for line in block:
                    write(line)
                count += 1
                as_min, as_max = min(as_min, alpha), max(as_max, alpha)
        if not closed or count == 0:
            raise ValueError('Missing closing LHE tag or no events')
        temp.replace(destination)
    except BaseException:
        temp.unlink(missing_ok=True)
        raise
    summary = dict(pdf='NNPDF31_lo_as_0130', member=0, lhapdf_id=315200,
                   scale='SCALUP = sqrt(sum(final-parton pT^2)/6), iqopt=1, qfac=1',
                   events=count, replaced=replaced, alpha_s_min=as_min, alpha_s_max=as_max,
                   original_sha256=original_hash.hexdigest(), corrected_sha256=corrected_hash.hexdigest())
    destination.with_suffix('.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--lhapdf-base', type=Path, required=True)
    args = parser.parse_args()
    convert(args.input, args.output, args.lhapdf_base)
