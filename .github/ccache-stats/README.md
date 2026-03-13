#  ccache-states

This action prints "ccache -s" to stdout inside the gh actions vms

## Inputs

N/A

## Outputs

### `stats`

The stats reported in `ccache -s` call.

## Build

Ensure you have `rollup` installed and available on the PATH. Then do
```bash
rollup --config rollup.config.js  --bundleConfigAsCjs
# or you can do
npm run build
```
This will generate the files in the expected location for gh action.
After these changes make a commit and open a PR.
