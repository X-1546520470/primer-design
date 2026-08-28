#!/bin/zsh
set -euo pipefail

APP_DIR="${0:A:h}"
APP_PYTHON="/opt/anaconda3/envs/Probe/bin/python"

if [[ ! -x "$APP_PYTHON" ]]; then
  print -u2 "找不到 Probe 环境：$APP_PYTHON"
  print -u2 "请先确认 /opt/anaconda3/envs/Probe 已安装。"
  read -r "?按回车键关闭…"
  exit 1
fi

export PYTHONNOUSERSITE=1
export PATH="/opt/anaconda3/envs/Probe/bin:$PATH"
cd "$APP_DIR"
exec "$APP_PYTHON" -m streamlit run "$APP_DIR/app.py"
