import os
import re
import shutil

# ================= USER CONFIG =================
BASE_DIR = "/Users/davidimig/projects/squirel/data/chains"

cases  = [
    'LCDM', 
    'DR', 
    'FD', 
    'BE', 
    'LN', 
    'LN_FD_fid', 
    'om_L_FD', 
    'fixed_zNR_FD', 
    'Tup_FD',
    'logzNR4_FD',
    'zup_FD',
    ]
cases = [cases[-1]]
probes = ['PR3', 'PR3+lens', 'PR3+eBOSS', 'PR4', 'S4', 'priors']
# probes = [probes[2]]
nus    = ['N_mnu=1_Mmnu=0.06', 'N_mnu=2_Mmnu=0.11', 'N_mnu=3_Mmnu=0.06']

# ===============================================


def process_txt_file(path):
    with open(path) as f:
        lines = f.readlines()

    n_lines = len(lines)

    update_idx = [
        i for i, l in enumerate(lines)
        if "update" in l
    ]

    if not update_idx:
        ans = input(
            f"No 'update' found in {os.path.basename(path)}"
            f" ({n_lines} total lines). "
            "Delete file? [y/N]: "
        ).strip().lower()

        if ans == "y":
            os.remove(path)
            print(f"    deleted {os.path.basename(path)} ({n_lines} lines)")
        else:
            print(f"    kept {os.path.basename(path)} unchanged")
        return

    cut = update_idx[-1]
    n_trimmed = cut

    with open(path, "w") as f:
        f.writelines(lines[cut:])

    print(
        f"    trimmed {os.path.basename(path)} "
        f"(removed {n_trimmed} lines, kept {n_lines - n_trimmed})"
    )


for case in cases:
    for probe in probes:
        for nu in nus:

            string2 = f"{case}_{probe}_{nu}"

            # if string2 != 'FD_PR3_N_mnu=1_Mmnu=0.06':
            #     continue

            directory = os.path.join(BASE_DIR, probe, nu, string2)

            if not os.path.isdir(directory):
                continue

            print(f"\nEntering directory: {directory}")

            # ---- rename ALL .txt files ----
            txt_files = sorted(
                n for n in os.listdir(directory)
                if n.endswith(".txt") and n != "log.param"
            )

            for i, old in enumerate(txt_files):
                new = f"{string2}__{i}.txt"
                os.rename(
                    os.path.join(directory, old),
                    os.path.join(directory, new),
                )
                print(f"  renamed {old} → {new}")

            # ---- handle paramnames ----
            param_files = [
                n for n in os.listdir(directory)
                if n.endswith(".paramnames")
            ]

            if param_files:
                keep = param_files[0]
                new_param = f"{string2}_.paramnames"
                os.rename(
                    os.path.join(directory, keep),
                    os.path.join(directory, new_param),
                )
                print(f"  kept paramnames → {new_param}")

                for extra in param_files[1:]:
                    os.remove(os.path.join(directory, extra))
                    print(f"  deleted extra paramnames {extra}")

            # ---- delete everything else ----
            for name in os.listdir(directory):
                path = os.path.join(directory, name)

                if os.path.isdir(path):
                    shutil.rmtree(path)
                    print(f"  removed directory {name}")
                    continue

                if (
                    name == "log.param"
                    or re.match(rf"^{re.escape(string2)}__\d+\.txt$", name)
                    or name == f"{string2}_.paramnames"
                ):
                    continue

                os.remove(path)
                print(f"  deleted {name}")

            # ---- ALWAYS trim contents of txt files ----
            for name in os.listdir(directory):
                if re.match(rf"^{re.escape(string2)}__\d+\.txt$", name):
                    process_txt_file(os.path.join(directory, name))
