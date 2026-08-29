"""MycoPrimerV2 启动入口：`python -m mycoprimer` 或控制台命令 `mycoprimer`。

数据目录解析顺序：
    1. 环境变量 PROBESTUDIO_HOME 指定的目录；
    2. 当前工作目录（若其中已有 genome_registry.json，即项目目录/既有数据目录）；
    3. ~/ProbeStudioData（首次在其它位置启动时自动创建）。

主题配置随包分发（config_default.toml），首次启动时复制到数据目录，
保证任何位置启动都有统一界面主题。

额外参数原样透传给 streamlit，例如：
    mycoprimer --server.port 8600
    mycoprimer --server.headless true
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def main() -> int:
    env_home = os.environ.get("PROBESTUDIO_HOME")
    if env_home:
        data_home = Path(env_home)
    elif (Path.cwd() / "genome_registry.json").is_file():
        data_home = Path.cwd()  # 项目目录或既有数据目录，保持原位
    else:
        data_home = Path.home() / "ProbeStudioData"
    data_home.mkdir(parents=True, exist_ok=True)

    # 主题自举：数据目录没有主题配置时，使用随包分发的默认主题。
    theme_dir = data_home / ".streamlit"
    theme_file = theme_dir / "config.toml"
    bundled = Path(__file__).with_name("config_default.toml")
    if bundled.is_file() and not theme_file.is_file():
        theme_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(bundled, theme_file)

    app_path = Path(__file__).with_name("app.py")
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print(
            "缺少 streamlit，请先安装依赖：pip install -e . （在项目目录内）",
            file=sys.stderr,
        )
        return 1

    sys.argv = ["streamlit", "run", str(app_path), *sys.argv[1:]]
    os.chdir(data_home)
    return stcli.main()


if __name__ == "__main__":
    raise SystemExit(main())
