#!/usr/bin/env python3
"""Export additive GenFilterInfo counters, including jobs with no accepted events."""
import argparse
import json
from pathlib import Path


def summarise(filename, mode):
    from DataFormats.FWLite import Handle, Lumis

    handle = Handle("GenFilterInfo")
    getters = {
        "n_before_positive": "numTotalPositiveEvents",
        "n_before_nonpositive": "numTotalNegativeEvents",
        "n_pass_positive": "numPassPositiveEvents",
        "n_pass_nonpositive": "numPassNegativeEvents",
        "sumw_before": "sumWeights",
        "sumw2_before": "sumWeights2",
        "sumw_pass": "sumPassWeights",
        "sumw2_pass": "sumPassWeights2",
    }
    result = {key: 0 for key in getters}
    n_lumis = 0
    for lumi in Lumis(filename):
        lumi.getByLabel("centralJetFilterEfficiency", handle)
        if not handle.isValid():
            raise RuntimeError("Missing centralJetFilterEfficiency lumi product")
        product = handle.product()
        for key, getter in getters.items():
            result[key] += getattr(product, getter)()
        n_lumis += 1
    if not n_lumis:
        raise RuntimeError("No luminosity blocks: cannot determine filter efficiency")
    result["n_before"] = result["n_before_positive"] + result["n_before_nonpositive"]
    result["n_pass"] = result["n_pass_positive"] + result["n_pass_nonpositive"]
    result["count_efficiency"] = (
        result["n_pass"] / result["n_before"] if result["n_before"] else None
    )
    result["weighted_efficiency"] = (
        result["sumw_pass"] / result["sumw_before"] if result["sumw_before"] else None
    )
    result.update(
        input=str(filename), mode=mode, luminosity_blocks=n_lumis,
        selection="at least 6 anti-kT R=0.4 visible gen jets, pt>20 GeV, abs(eta)<3.0",
        denominator="events with generator GenEventInfoProduct after shower/merging",
        note="off records the would-pass efficiency without rejecting events; sum counters across jobs before dividing",
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=["on", "off"], required=True)
    args = parser.parse_args()
    text = json.dumps(summarise(args.input, args.mode), indent=2, allow_nan=False) + "\n"
    args.output.write_text(text)
    print(text, end="")
