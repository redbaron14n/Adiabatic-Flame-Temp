# Adiabatic-Flame-Temp

## UV package manager

Install UV on your PC.

```
pip install uv
```

In project directory run

```
uv sync
```

to create a populate .venv folder with the libraries used by the project.
To activate the virtual environment, run

```
source .venv/bin/activate
```

To deactivate the virtual environment, run

```
deactivate
```

To select the venv as the workspace venv, use ctl-shift-p 'python: select interpreter' and select the .venv/scripts/python.exe.

To run the app

```
uv run python -m app
```
