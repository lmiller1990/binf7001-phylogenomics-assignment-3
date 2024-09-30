import os
from tabulate import tabulate

dirs = os.listdir("contig_stats")


def get_percent_assembly(kmer):
    with open(os.path.join("velvet", kmer, "Log"), "r") as f:
        for line in f.readlines():
            if "Final graph" in line:
                _, content = line.split("using")
                total, rest = content.split("reads")
                num, denom = map(lambda x: int(x), total.split("/"))
                return f"{num} / {denom} ({round((num/denom * 100), 2)}%)"


def parse_contig_stats(path):
    stats = []

    with open(os.path.join("contig_stats", path), "r") as f:
        for line in f.readlines():
            if "Main genome scaffold N/L50" in line:
                stats.append(["Scaffold N/L50", line.split(":")[-1].strip()])

            if "Max scaffold len" in line:
                stats.append(["Max scaffold length", line.split(":")[-1].strip()])

            if "# of scaffolds > 50 KB" in line:
                stats.append(["Number of scaffolds > 50 KB", line.split(":")[-1].strip()])

            if "All" in line:
                n = (line.split("\t")[3].strip())
                stats.append(["Total scaffold length", n])

    s = path.split('.contig')
    kmer = s[0].strip()
    total = get_percent_assembly(kmer)
    stats.append(["% assembled", total])

    return stats

def collect():
    kmer_data = []
    for dir in sorted(dirs):
        stats = parse_contig_stats(dir)
        stats.insert(0, ["kmer", dir.split(".")[0]])
        kmer_data.append(stats)
    return kmer_data

if __name__ == "__main__":
    data = collect()
    rows = list(map(lambda x: x[0], data[0]))
    header = "| " + " | ".join(rows) + " |"
    header += "\n"
    header += "|" + "".join(["------|" for x in range(len(rows))])

    lines = []
    for row in data:
        line = "| "
        for item in row:
            line += item[1] + " | "
        lines.append(line)

    all = "\n".join(lines)
    print(header + "\n" + all)


