# Agent Instructions

## Code Style

- Use the fewest lines that still satisfy correctness, APIs, framework contracts, and readability.
- Do not add filler comments, narrative docstrings on obvious symbols, redundant guards, or try/except blocks that do not handle anything meaningful.
- Prefer early returns, comprehensions or generator expressions, inline literals, and single-expression helpers when readability remains clear.
- When editing, touch only what is needed. Do not wrap simple logic in extra layers just to explain a change.

## Python Scripts

- Do not add command-line interfaces or `argparse` to new Python scripts unless explicitly requested.
- Prefer editable module-level constants near the top of the script for run options such as paths, grids, modes, and output filenames.
- Module-level constants/global run settings are acceptable and preferred for generated scripts; do not avoid them when they make local reruns easier.
- Do not hard-code aircraft architecture counts such as number of engines, propulsive motors, generators, turbines, motors per engine, or generators per motor in scripts. Store those as aircraft/propulsion object attributes and pass the aircraft or propulsion object into helpers.
- Keep scripts directly runnable with `python path/to/script.py`.
- If a script needs reusable behavior, expose functions and keep `main()` as a thin call using the constants.
- Existing CLI scripts may keep their current interface, but do not expand CLI arguments without explicit approval.

## Data Flow And State

- Pass data explicitly through function arguments, return values, constructors, small context objects, or framework inputs/outputs.
- Module-level constants are fine when immutable by convention, especially script run options. Avoid runtime mutation of module-level "constants".
- Do not use globals for cross-cutting mutable state such as caches, current config, last result, or debug flags. Use parameters, object attributes, dedicated instances, or framework hooks.
- Prefer class or instance attributes when state logically belongs to an object. Prefer dependency injection over importing a mutable singleton.
- When integrating with OpenMDAO-style models, wire values through the model graph or problem setup instead of module globals that components read implicitly.
