"""ProbeStudio 启动入口：`python -m probestudio` 或控制台命令 `probestudio`。

行为：
    1. 确定数据目录：环境变量 PROBESTUDIO_HOME，默认 ~/ProbeStudioData；
       双击 launch.command 启动时其工作目录已是项目目录，数据仍在项目内。
    2. 在该目录下启动 Streamlit 服务（本地 127.0.0.1，序列不出机）。

额外参数原样透传给 streamlit，例如：
    probestudio --server.port 8600
    probestudio --server.headless true
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> int:
    default_home = Path.home() / "ProbeStudioData"
    data_home = Path(os.environ.get("PROBESTUDIO_HOME", default_home))
    data_home.mkdir(parents=True, exist_ok=True)
    os.chdir(data_home)

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
    return stcli.main()


if __name__ == "__main__":
    raise SystemExit(main())
