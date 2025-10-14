import re
from datetime import datetime
import pytz
import pandas as pd
import matplotlib.pyplot as plt
from dateutil import parser


py_log_file = True

if py_log_file:
    log_dir = "Z:/decentralized_choice/code/2_structural_estimation/logs/estimation_15115451.log"

    cuda = False

    
    with open(log_dir, "r") as f:
        content = f.read()

    # Count number of times objfn values were calculated
    count_objfn = content.count("Objfn Value")
    count_iter = content.count("At iterate")


    # Extract the start time using regex
    match = re.search(r"Starting estimation at (.+)", content)
    if match:
        start_time_str = match.group(1).strip()

        # Remove the timezone part (e.g., 'CDT') and parse the rest
        parts = start_time_str.split()
        cleaned_str = ' '.join(parts[:4] + parts[5:])  # drop timezone (index 4)
        start_time_naive = datetime.strptime(cleaned_str, "%a %b %d %H:%M:%S %Y")

        # Attach Central Timezone (CDT)
        central = pytz.timezone("US/Central")
        start_time = central.localize(start_time_naive)

       # Get current time
        now = datetime.now(central)

        # Extract the time string
        match = re.search(r"Job completed at (.+)", content)

        if match:
            date_str = match.group(1).strip()
            end_time = parser.parse(date_str)  # Automatically handles spaces, non-zero-padded day, etc.
            # Make sure it's timezone-aware (often it already is, but double-check)
            if end_time.tzinfo is None:
                end_time = central.localize(end_time)
            last_time = end_time
        else:
            last_time = now

        # Compute elapsed seconds
        elapsed_seconds = (last_time - start_time).total_seconds()

        # Compute iteration rate
        rate = count_objfn / elapsed_seconds *60 if elapsed_seconds > 0 else float('inf')

        print(f"Estimation started at: {start_time}")
        print(f"Last time: {last_time}")
        print(f"Elapsed time: {elapsed_seconds:.2f} seconds")
        print(f"Objfn evaluations: {count_objfn}")
        print(f"Rate: {rate:.4f} evaluations per minute")
        print(f"Iterations: {count_iter}")
        print(f"Rate: {count_iter/elapsed_seconds *60}")
    else:
        print("Start time not found in the log.")


    if not cuda:
        matches = re.findall(r"At iterate\s+(\d+)\s+f=\s+([0-9D\+\.\-]+)", content)
    else:
        matches = re.findall(r"Objfn Value:\s*([0-9.eE+-]+)", content)


    print(matches)
    if not cuda:
        df = pd.DataFrame(matches, columns=["iter", "f"])
        #df  = df.iloc[200:]

    else:
        df=  pd.DataFrame(matches, columns=["f"])
        df["iter"] = df.index
        df  = df.loc[200:,:]
  

   # Convert to correct types
    df["iter"] = df["iter"].astype(int)
    df["f"] = df["f"].astype(str).str.replace("D", "E").astype(float)
   # df = df[df["iter"]>=2000]

    # Identify when iteration resets
    df["attempt"] = (df["iter"].diff() < 0).cumsum() + 1  # Start from attempt 1

    # Plot each attempt as a separate line
    plt.figure(figsize=(10, 5))
    for attempt, subset in df.groupby("attempt"):
        plt.plot(subset["iter"], subset["f"], marker='o', label=f"Attempt {attempt}")

    plt.xlabel("Iteration")
    plt.ylabel("Objective Function Value (f)")
    plt.title("Objective Function over Iterations")
    plt.grid(True)
    plt.legend()
    plt.show()

    if not df.empty:
        f_start = df["f"].iloc[0]
        f_end = df["f"].iloc[-1]
        f_drop = f_start - f_end
        print(f"\nObjective function decreased by: {f_drop:.4f}")
        print(f"Average decrease per minute: {f_drop / elapsed_seconds * 60:.4f}")

else:
    log_dir = "Z:/decentralized_choice/eta_out_K_3_start2025_07_31_12_26_14.txt"

    with open(log_dir, "r") as f:
        content = f.read()

    # Extract datetime and value pairs
    matches = re.findall(r"\[objfn\s*-\s*(.*?)\]\s*Objective function value:\s*([0-9eE\+\.\-]+)", content)

    # Build DataFrame
    df_objfn_time = pd.DataFrame(matches, columns=["timestamp", "objfn_value"])
    df_objfn_time["timestamp"] = pd.to_datetime(df_objfn_time["timestamp"])
    df_objfn_time["objfn_value"] = df_objfn_time["objfn_value"].astype(float)
    df_objfn_time = df_objfn_time.loc[10:,]
    print(df_objfn_time)


    # Plot original and 20-MA
    plt.figure(figsize=(10, 5))
    plt.plot(df_objfn_time["timestamp"], df_objfn_time["objfn_value"], label="Objective Function", marker='o')
    # plt.plot(df_objfn_time["timestamp"], df_objfn_time["objfn_20ma"], label="20-MA", linestyle='--', color='orange')

    plt.xlabel("Timestamp")
    plt.ylabel("Objective Function Value")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()