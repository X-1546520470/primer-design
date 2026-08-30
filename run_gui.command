#!/bin/zsh
# MycoPrimerV2 桌面版双击启动脚本（Tkinter 原生窗口，无需浏览器）。
# 部署：首次使用请先双击 setup.command。
set -euo pipefail
APP_DIR="${0:A:h}"
cd "$APP_DIR"
source "$APP_DIR/_find_python.sh"
exec "$APP_PYTHON" mycoprimer_gui.py
