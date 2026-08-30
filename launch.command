#!/bin/zsh
# MycoPrimerV2 网页版双击启动脚本（Streamlit，功能与桌面版相同）。
# 部署：首次使用请先双击 setup.command。
set -euo pipefail
APP_DIR="${0:A:h}"
cd "$APP_DIR"
source "$APP_DIR/_find_python.sh"
exec "$APP_PYTHON" -m mycoprimer
