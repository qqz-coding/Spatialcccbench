import re
import os

def find_time_1(path,line):

    # 检查文件是否存在
    if not os.path.isfile(path):
        print(f"错误：文件不存在 - {path}")
        return "NA"

    # 检查文件是否可读
    if not os.access(path, os.R_OK):
        print(f"错误：没有读取权限 - {path}")
        return "NA"

    # 检查文件是否为空
    if os.path.getsize(path) == 0:
        print(f"警告：文件为空 - {path}")
        return "NA"

    try:
        with open(path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

            if not lines:
                print("警告：文件内容为空")
                return "NA"

            last_line = lines[-2].strip()

            if not last_line:
                print("警告：最后一行为空")
                return "NA"

            # 更健壮的正则表达式，匹配各种数字格式
            # 匹配: 整数、浮点数、科学计数法
            numbers_found = re.findall(r'-?\d*\.?\d+(?:[eE][-+]?\d+)?', last_line)

            if numbers_found:
                try:
                    # 取第一个匹配的数字并转换为float
                    time_value = float(numbers_found[0])
                    print(f"成功提取时间: {time_value} (从: '{last_line}')")
                    return time_value
                except (ValueError, TypeError) as e:
                    print(f"警告：数字格式错误 '{numbers_found[0]}' - {e}")
                    return "NA"
            else:
                print(f"警告：未找到数字格式内容。最后一行: '{last_line}'")
                return "NA"

    except FileNotFoundError:
        print(f"错误：文件不存在 - {path}")
        return "NA"
    except PermissionError:
        print(f"错误：没有文件读取权限 - {path}")
        return "NA"
    except UnicodeDecodeError:
        # 尝试其他编码
        try:
            with open(path, 'r', encoding='latin-1') as file:
                lines = file.readlines()
                last_line = lines[-2].strip() if lines else ""
                numbers_found = re.findall(r'-?\d*\.?\d+(?:[eE][-+]?\d+)?', last_line)
                if numbers_found:
                    return float(numbers_found[0])
                return "NA"
        except Exception as e:
            print(f"错误：文件编码问题 - {e}")
            return "NA"
    except Exception as e:
        print(f"错误：读取文件时发生未知错误 - {e}")
        return "NA"
