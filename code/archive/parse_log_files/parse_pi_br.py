
import os, re
import pandas as pd
from itertools import chain

log_path = os.path.normpath ("Z:/decentralized_choice/code/logs")
out_path = os.path.normpath("Z:/decentralized_choice/estimates")


float_pat = r"""
[-+]?                # optional sign
(?:\d+\.\d+|\d+)     # number with or without decimal
(?:[eE][-+]?\d+)?    # optional exponent
"""


def extract_model_k_info(text: str) -> dict:

    """
    Look for lines like:
      [1] "Model: eta_out"
      [1] "K: 3"
    and return {"model": "eta_out", "K": 3}.
    """
    # strip R's [1] prefix and quotes
    cleaned = re.sub(r'\[\d+\]\s*"', '', text)
    cleaned = cleaned.replace('"', '')

    # find "Model: ..." and "K: ..."
    model_match = re.search(r"Model:\s*(\S+)", cleaned)
    k_match = re.search(r"K:\s*(\d+)", cleaned)

    info = {}
    if model_match:
        info["model"] = model_match.group(1)
    if k_match:
        info["K"] = int(k_match.group(1))


    scenarios = [
    "Shutting Down Application Costs",
    "Info Provision",
    "Targeted Info Provision"
    ]
    for s in scenarios:
         if re.search(s, cleaned):
             info["scenario"] = s

    return info

def parse_blocks_backward(log_text: str) -> pd.DataFrame:
    CLOSE_LINE = re.compile(r'^\s*\[\s*396\]')
    OPEN_LINE  = re.compile(r'^\s*\[\s*1\]')

    lines = log_text.splitlines()
    i = len(lines) - 1

    nums_ls = []

    while i >= 0:
        if CLOSE_LINE.match(lines[i]):
            j = i
            while j >= 0 and not OPEN_LINE.match(lines[j]):
                j -= 1

            if j >= 0:
                def parse_numbers(s: str):
                    # remove [NNN] safely with spaces so tokens don't glue
                    s = re.sub(r'\[\s*\d+\]', ' ', s)

                    STRICT_NUM = re.compile(r'[+-]?\d\.\d{6}e[+-]\d{2}')
                    return tuple(float(x) for x in re.findall(STRICT_NUM, s))

                per_line = [parse_numbers(lines[k]) for k in range(j, i + 1)]
                flat = tuple(chain.from_iterable(per_line))  # <- flatten tuples

                nums_ls.insert(0, flat)  # keep natural (top-to-bottom) order
                
                #if max(["Gradient Norm:[296]" in l for l in lines[j:i+1]]):
                #if j == 134064:
                #    print(flat[310:320])
                #    return lines[j:i+1]
                
                i = j - 1
                continue
        i -= 1

    # remove duplicate blocks by converting to dict (preserves order)
    unique_blocks = list(dict.fromkeys(nums_ls))

    # expand into DataFrame
    rows = []
    for block_id, block in enumerate(unique_blocks, start=1):
        rows.extend(
            {"pi_br": val, "pi_idx": pos, "block_idx": block_id}
            for pos, val in enumerate(block, start=1)
        )

    return pd.DataFrame(rows)


log_path = os.path.normpath("Z:/decentralized_choice/code/logs")
log_ls = [15000472, 15000473, 15000474]  # or your full list

all_results = []
for l in log_ls:
    f_name = os.path.join(log_path, f"counterfactual_{l}.log")
    with open(f_name, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
        text = text.replace("Current Gradient Norm:", "")
        text = text.replace("Round", "")

        info = extract_model_k_info(text)
        results = parse_blocks_backward(text)
        results["K"] = info["K"]
        results["model"] = info["model"]
        results["scenario"] = info["scenario"]
        results["logfile"] = l
        all_results.append(results)


combined = pd.concat(all_results, ignore_index=True)
combined.to_csv(os.path.join(out_path, "pi_br_by_log_and_scenario.csv"), index=False)
print(combined)

# Group by (scenario, model, K, pi_idx)
summary = (
    combined.groupby(["scenario", "model", "K", "pi_idx"])
            .agg(
                mean_pi_br=("pi_br", "mean"),
                n_unique_logs=("block_idx", "nunique")
            )
            .reset_index()
)

# Optional: look
print(summary.head())

summary.to_csv(os.path.join(out_path,"pi_br_average.csv"), index=False)
