#!/bin/zsh
# setup.command — MycoPrimerV2 一键部署（在新电脑上 clone/下载后先双击这个）。
#
# 部署策略（两条自洽路径，二选一）：
#   A. 检测到 conda：创建专用环境 mycoprimer（python 3.11 + bowtie2），
#      并在该环境内 pip install 本项目 —— tkinter/tk 随 conda 自带，一站式。
#   B. 无 conda 但有 Homebrew：brew 安装 python@3.12 + python-tk + bowtie2，
#      再用项目内 .venv 安装 Python 依赖。
# 都没有：打印手动安装指引。
#
# 部署成功后把解释器路径写入 .python_path，启动脚本据此运行。
# 完成后双击 run_gui.command（桌面版）或 launch.command（网页版）。

set -euo pipefail

APP_DIR="${0:A:h}"
cd "$APP_DIR"

print "==> MycoPrimerV2 部署开始（目录：$APP_DIR）"

_ok() { "$1" -c "import primer3, Bio" >/dev/null 2>&1; }

# ------------------------------------------------------------------
# 路径 A：conda（推荐——一个环境装齐 python + tkinter + bowtie2 + 依赖）
# ------------------------------------------------------------------
if command -v conda >/dev/null 2>&1; then
  print "==> 检测到 conda，创建专用环境 mycoprimer（python 3.11 + bowtie2）…"
  conda create -n mycoprimer -c bioconda -c conda-forge \
      python=3.11 bowtie2 tk -y
  conda_base="$(conda info --base)"
  env_python="$conda_base/envs/mycoprimer/bin/python"
  print "==> 在 mycoprimer 环境内安装 MycoPrimerV2 及 Python 依赖…"
  "$env_python" -m pip install --upgrade pip -q
  "$env_python" -m pip install -e "$APP_DIR" -q
  print "$env_python" > .python_path
  if _ok "$env_python"; then
    print ""
    print "✅ 部署完成（conda 路径）。双击 run_gui.command 即可使用。"
    read -r "?按回车键关闭…"
    exit 0
  fi
  print "⚠️  conda 路径校验未通过，尝试 Homebrew 路径…"
fi

# ------------------------------------------------------------------
# 路径 B：Homebrew（python + python-tk + bowtie2）+ 项目内 .venv
# ------------------------------------------------------------------
if command -v brew >/dev/null 2>&1; then
  print "==> 使用 Homebrew 安装 python@3.12 / python-tk / bowtie2…"
  brew install --quiet python@3.12 python-tk@3.12 bowtie2 || true
  brew_bin="$(brew --prefix)/bin"
  export PATH="$brew_bin:$PATH"
fi

# 选一个 ≥3.10 的解释器创建 .venv
PY=""
for cand in python3.12 python3.11 python3.10 python3 \
            /opt/homebrew/bin/python3.12 /usr/local/bin/python3.12 /usr/bin/python3; do
  if command -v "$cand" >/dev/null 2>&1; then
    ver="$("$cand" -c 'import sys; print("%d.%d" % sys.version_info[:2])' 2>/dev/null || echo 0)"
    if [[ "$($cand -c 'import sys; print(sys.version_info[0])')" == "3" ]] \
       && [[ "${ver#*.}" -ge 10 ]]; then
      PY="$(command -v "$cand")"
      print "==> 找到 Python $ver：$PY"
      break
    fi
  fi
done
if [[ -z "$PY" ]]; then
  print -u2 "✗ 未找到 Python ≥ 3.10。请安装 Miniconda 或 Homebrew Python 后重试。"
  read -r "?按回车键关闭…"
  exit 1
fi

if [[ ! -x ".venv/bin/python" ]]; then
  print "==> 创建项目内虚拟环境 .venv…"
  "$PY" -m venv .venv
fi
print "==> 安装 Python 依赖（可能需要几分钟）…"
.venv/bin/python -m pip install --upgrade pip -q
.venv/bin/python -m pip install -e . -q
print "$APP_DIR/.venv/bin/python" > .python_path

if _ok "$APP_DIR/.venv/bin/python"; then
  if command -v bowtie2 >/dev/null 2>&1 || [[ -x "$brew_bin/bowtie2" ]] \
     || [[ -x "/opt/homebrew/bin/bowtie2" ]]; then
    print ""
    print "✅ 部署完成（Homebrew 路径）。双击 run_gui.command 即可使用。"
  else
    print ""
    print "✅ Python 依赖完成；⚠️ bowtie2 尚未就绪——桌面版可打开，"
    print "   但涉及比对的过滤会报错。请执行：brew install bowtie2"
  fi
else
  print -u2 "⚠️ 依赖校验未通过，请把上方报错发给维护者。"
fi
read -r "?按回车键关闭…"
