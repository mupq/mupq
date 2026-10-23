# mupq
Provides common files for instances of **mupq**, e.g., for 
 - https://github.com/mupq/pqm4
 - https://github.com/mupq/pqm3 
 - https://github.com/mupq/pqriscv

# Several boards on one host
When more than one debug adapter is attached, the flash step has to be told which one to use. Two environment variables select the adapter by its serial number:

| Variable | Used by | Effect |
|---|---|---|
| `MUPQ_OPENOCD_SERIAL` | `OpenOCD`, `OpenOCDITM` | passed to OpenOCD as `adapter serial <serial>` |
| `MUPQ_ST_FLASH_ARGS` | `StLink` | extra arguments to `st-flash`, e.g. `--serial <serial>` |

Unset, the flash commands are unchanged.

The serial numbers of the attached ST-LINK adapters are listed by
```
st-info --probe
```
or, on macOS, by
```
system_profiler SPUSBHostDataType
```

For example, to flash a board through OpenOCD and another through st-flash on the same host:
```
export MUPQ_OPENOCD_SERIAL=066BFF575550897767144107
export MUPQ_ST_FLASH_ARGS="--serial 066CFF485648756687012533"
```

# License
Different parts of **mupq** have different licenses. Each subdirectory containing implementations contains a LICENSE file stating under what license that specific implementation is released. The files in common contain licensing information at the top of the file (and are currently either public domain or MIT).
All other code in this repository is dual-licensed under [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0) and under the conditions of [CC0](https://creativecommons.org/publicdomain/zero/1.0/).
