#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--itf", type=str, required=True, help="Path to the raw ITF file", metavar="PATH")
    parser.add_argument("--skeleton", type=str, required=True, help="Path to the skeleton file", metavar="PATH")

    args = parser.parse_args()

    raw = open(args.itf, "r").read()
    skeleton = open(args.skeleton, "r").read()

    # Tensor renaming
    raw = raw.replace("HAM_D:cc", "g:cc")
    raw = raw.replace("HAM_D:ee", "g:ee")
    raw = raw.replace("HAM_D:aaaa", "K:aaaa")
    raw = raw.replace("HAM_D:aa", "f:aa")
    raw = raw.replace("GAM0:aaaa", "Ym2")
    raw = raw.replace("GAM0:aa", "Ym1")
    raw = raw.replace("T2g:", "T2:")

    # Strip ---- end (and anything below it)
    raw = raw[:raw.index("---- end")].strip()
    skeleton = skeleton[:skeleton.index("---- end")].strip()

    declarations = skeleton[:skeleton.index("---- code")].strip()
    raw_declarations = raw[raw.index("---- decl") : raw.index("---- code")].strip()

    for line in raw_declarations.splitlines():
        if not line.startswith("tensor:"):
            continue

        tensor_spec = line[:line.index("[") + 1]
        
        if not tensor_spec in declarations:
            declarations += "\n" + line

    itf = declarations
    itf += "\n"
    itf += skeleton[skeleton.index("---- code"):]
    itf += "\n\n"
    itf += raw[raw.index("---- code"):]
    itf += "\n"

    if not "code(\"Update_INTkx\")" in itf:
        itf += "\n---- code(\"Update_INTkx\")\n"
        itf += "# Set INTkx tensors to zero\n"
        itf += "alloc INTkx:eeaa[abuv]\nstore INTkx:eeaa[abuv]\n\n"
        itf += "alloc INTkx:eeac[abui]\nstore INTkx:eeac[abui]\n\n"
        itf += "alloc INTkx:eecc[abij]\nstore INTkx:eecc[abij]\n\n"

    if not "code(\"Sum_T1\")" in itf:
        itf += "\n---- code(\"Sum_T1\")\n"
        itf += "alloc T1s:ec[ai]\n"
        itf += "load T1:ec[ai]\n"
        itf += ".T1s:ec[ai] += T1:ec[ai]\n"
        itf += "drop T1:ec[ai]\n"
        itf += "load T2:ec[ai]\n"
        itf += ".T1s:ec[ai] += T2:ec[ai]\n"
        itf += "drop T2:ec[ai]\n"
        itf += "store T1s:ec[ai]\n"

    itf += "\n\n---- end\n"

    print(itf)


if __name__ == "__main__":
    main()

