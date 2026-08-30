#!/bin/zsh
# _find_python.sh — 解析 MycoPrimerV2 使用的 Python 解释器（被各启动脚本 source）。
#
# 解析顺序（找到第一个"能 import 引擎依赖"的解释器即返回）：
#   1. 环境变量 PROBESTUDIO_PYTHON
#   2. 项目内虚拟环境 .venv/bin/python
#   3. 记录文件 .python_path（setup.command 部署成功后写入）
#   4. 本机既有 conda 环境 Probe（历史路径，向前兼容）
#   5. PATH 中的 python3
#
# 用法：source 本文件后，使用变量 APP_PYTHON；失败时脚本已打印指引并 exit 1。

APP_DIR="${0:A:h}"

_python_ok() {
  # 能否 import 引擎的核心依赖（primer3 + Bio）
  "$1" -c "import primer3, Bio" >/dev/null 2>&1
}

_resolve() {
  # 1. 环境变量
  if [[ -n "${PROBESTUDIO_PYTHON:-}" && -x "$PROBESTUDIO_PYTHON" ]] && _python_ok "$PROBESTUDIO_PYTHON"; then
    print -r -- "$PROBESTUDIO_PYTHON"
    return 0
  fi
  # 2. 项目内虚拟环境
  if [[ -x "$APP_DIR/.venv/bin/python" ]] && _python_ok "$APP_DIR/.venv/bin/python"; then
    print -r -- "$APP_DIR/.venv/bin/python"
    return 0
  fi
  # 3. 部署记录文件
  if [[ -f "$APP_DIR/.python_path" ]]; then
    local recorded
    recorded="$(cat "$APP_DIR/.python_path" 2>/dev/null | head -1)"
    if [[ -x "$recorded" ]] && _python_ok "$recorded"; then
      print -r -- "$recorded"
      return 0
    fi
  fi
  # 4. 历史路径（本机 conda 环境）
  if [[ -x "/opt/anaconda3/envs/Probe/bin/python" ]] && _python_ok "/opt/anaconda3/envs/Probe/bin/python"; then
    print -r -- "/opt/anaconda3/envs/Probe/bin/python"
    return 0
  fi
  # 5. PATH 中的 python3
  if command -v python3 >/dev/null 2>&1 && _python_ok "$(command -v python3)"; then
    print -r -- "$(command -v python3)"
    return 0
  fi
  return 1
}

if APP_PYTHON="$(_resolve)"; then
  export APP_PYTHON
  # bowtie2 可能不在 PATH：把解释器所在环境与常见位置补进 PATH
  export PATH="${APP_PYTHON:h}:$PATH"
else
  print -u2 "MycoPrimerV2 尚未部署：找不到可用的 Python（缺 primer3/biopython 依赖）。"
  print -u2 "请先双击 setup.command 完成一键部署，或手动执行："
  print -u2 "  python3 -m venv .venv && .venv/bin/pip install -e ."
  print -u2 "  并安装 bowtie2（brew install bowtie2 或 conda install -c bioconda bowtie2）"
  read -r "?按回车键关闭…"
  exit 1
fi
