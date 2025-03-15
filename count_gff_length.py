# 统计gff文件中，某一类型元素的总长度
import argparse

def calculate_total_length_and_count(gff_file, category):
    total_length = 0
    count = 0

    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#"):  # 跳过注释行
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # 确保是有效的 GFF 行
            feature_type = fields[2]  # 第三列是元素类型
            start, end = int(fields[3]), int(fields[4])  # 第四、五列是起止坐标
            
            if feature_type == category:
                count += 1
                total_length += end - start + 1  # 计算区间长度
    
    return count, total_length

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate total length and count of a specific feature in a GFF file.")
    parser.add_argument("-g", "--gff", required=True, help="Input GFF file")
    parser.add_argument("-c", "--category", required=True, help="Feature type to calculate length and count (e.g., gene, mRNA, miRNA)")
    
    args = parser.parse_args()
    count, total_length = calculate_total_length_and_count(args.gff, args.category)
    
    print(f"The total count of {args.category}: {count}")
    print(f"The total length of {args.category}: {total_length}")
