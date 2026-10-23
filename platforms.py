# SPDX-License-Identifier: Apache-2.0 or CC0-1.0
from mupq import mupq

import abc
import re
import serial
import subprocess
import tempfile
import time
import os
import tqdm

try:
    import chipwhisperer as cw
except ImportError:
    pass


class Qemu(mupq.Platform):

    start_pat = re.compile('.*={4,}\n', re.DOTALL)
    end_pat = re.compile('#\n', re.DOTALL)

    def __init__(self, qemu, machine):
        super().__init__()
        self.qemu = qemu
        self.machine = machine
        self.platformname = "qemu"

    def __enter__(self):
        return super().__enter__()

    def __exit__(self, *args, **kwargs):
        return super().__exit__(*args, **kwargs)

    def run(self, binary_path, expiterations=1):
        if expiterations > 1:
            pb = tqdm.tqdm(total=expiterations, leave=False, desc="Running...")
        args = [
            self.qemu,
            "-M",
            self.machine,
            "-nographic",
            "-semihosting",
            "-kernel",
            binary_path,
        ]
        self.log.info(f'Running QEMU: {" ".join(args)}')
        try:
            proc = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, encoding="ascii")
            output = ""
            while "#" not in output:
                buf = proc.stdout.readline()
                # print(buf)
                output += buf
                if expiterations > 1:
                    if "+" in buf:
                        pb.update(buf.count("+"))
                    else:
                        pb.refresh()
            proc.wait()
        except Exception:
            try:
                if proc and proc.poll is not None:
                    proc.kill()
            except Exception:
                pass
            raise
        start = self.start_pat.search(output)
        end = self.end_pat.search(output, start.end())
        if end is None:
            return 'ERROR'
        proc.wait()
        if expiterations > 1:
            pb.close()
        return output[start.end():end.start()]


class SerialCommsPlatform(mupq.Platform):

    # Start pattern is at least five equal signs
    start_pat = re.compile(b'.*={4,}\n', re.DOTALL)

    def __init__(self, tty="/dev/ttyACM0", baud=38400, timeout=1):
        super().__init__()
        self._dev = serial.Serial(tty, baud, timeout=timeout)

    def __enter__(self):
        return super().__enter__()

    def __exit__(self, *args, **kwargs):
        self._dev.close()
        return super().__exit__(*args, **kwargs)

    def run(self, binary_path, expiterations=1):
        if expiterations > 1:
            pb = tqdm.tqdm(total=expiterations, leave=False, desc="Running...")
        self._dev.reset_input_buffer()
        self.flash(binary_path)
        # Wait for the first equal sign
        if self._dev.read_until(b'=')[-1] != b'='[0]:
            raise RuntimeError('Timout waiting for start')
        # Wait for the end of the equal delimiter
        start = self._dev.read_until(b'\n')
        self.log.debug(f'Found start pattern: {start}')
        if self.start_pat.fullmatch(start) is None:
            raise RuntimeError('Start does not match')
        # Wait for the end
        output = bytearray()
        while len(output) == 0 or output[-1] != b'#'[0]:
            data = self._dev.read_until(b'#', 128)
            if expiterations > 1:
                if b"+" in data:
                    pb.update(data.count(b"+"))
                else:
                    pb.refresh()
            output.extend(data)
        if expiterations > 1:
            pb.close()
        return output[:-1].decode('utf-8', 'ignore')

    @abc.abstractmethod
    def flash(self, binary_path):
        pass


class OpenOCD(SerialCommsPlatform):
    """Flashes with OpenOCD. MUPQ_OPENOCD_SERIAL selects the adapter by its
    serial number."""

    def __init__(self, script, tty="/dev/ttyACM0", baud=38400, timeout=60):
        super().__init__(tty, baud, timeout)
        self.script = script

    def flash(self, binary_path):
        args = ["openocd", "-f", self.script]
        if os.getenv("MUPQ_OPENOCD_SERIAL") is not None:
            args += ["-c", f"adapter serial {os.getenv('MUPQ_OPENOCD_SERIAL')}"]
        args += ["-c", f"program {binary_path} verify reset exit"]
        subprocess.check_call(
            args,
            # stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


class StLink(SerialCommsPlatform):
    def flash(self, binary_path):
        extraargs = []
        if os.getenv("MUPQ_ST_FLASH_ARGS") is not None:
            extraargs = os.getenv("MUPQ_ST_FLASH_ARGS").split()
        subprocess.check_call(
            ["st-flash"] + extraargs + ["--reset", "write", binary_path, "0x8000000"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def itm_port0(data):
    """Decode ITM software packets, returning (port 0 payload, bytes consumed).

    A trailing incomplete packet is left unconsumed."""
    out = bytearray()
    i = 0
    while i < len(data):
        header = data[i]
        size = {1: 1, 2: 2, 3: 4}.get(header & 0x03)
        if size is None or header & 0x04:
            # Synchronization, overflow or hardware source packets: skip the header
            i += 1
            continue
        if i + 1 + size > len(data):
            break
        if header >> 3 == 0:
            out.extend(data[i + 1:i + 1 + size])
        i += 1 + size
    return bytes(out), i


class OpenOCDITM(mupq.Platform):
    """Flashes with OpenOCD and reads the output from ITM stimulus port 0
    through the SWO pin, captured by the same OpenOCD process.

    The firmware sets up the ITM and the SWO prescaler for its own clock;
    OpenOCD captures the pin at pin_freq. MUPQ_OPENOCD_SERIAL selects the
    adapter by its serial number."""

    # Start pattern is at least five equal signs
    start_pat = re.compile(b'.*={4,}\n', re.DOTALL)

    def __init__(self, script, tpiu, pin_freq=2000000, traceclk=24000000, timeout=60):
        super().__init__()
        self.script = script
        self.tpiu = tpiu
        self.pin_freq = pin_freq
        self.traceclk = traceclk
        self.timeout = timeout

    def run(self, binary_path, expiterations=1):
        if expiterations > 1:
            pb = tqdm.tqdm(total=expiterations, leave=False, desc="Running...")
        address = " 0x08000000" if binary_path.endswith(".bin") else ""
        with tempfile.TemporaryDirectory() as tmp:
            trace = os.path.join(tmp, "swo.bin")
            commands = [
                "init",
                f"program {binary_path}{address} verify",
                "reset halt",
                f"{self.tpiu} configure -protocol uart -output {trace} -traceclk {self.traceclk} -pin-freq {self.pin_freq} -formatter off",
                f"{self.tpiu} enable",
                "itm port 0 on",
                "resume",
            ]
            args = ["openocd", "-f", self.script]
            if os.getenv("MUPQ_OPENOCD_SERIAL") is not None:
                args += ["-c", f"adapter serial {os.getenv('MUPQ_OPENOCD_SERIAL')}"]
            for command in commands:
                args += ["-c", command]
            self.log.info(f'Running OpenOCD: {" ".join(args)}')
            proc = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            try:
                raw = bytearray()
                text = bytearray()
                start = None
                last = time.time()
                # OpenOCD creates the trace file when the TPIU is enabled
                while not os.path.exists(trace):
                    if proc.poll() is not None:
                        raise RuntimeError(f"OpenOCD exited with {proc.returncode}")
                    if time.time() - last > self.timeout:
                        raise RuntimeError("Timeout waiting for the trace file")
                    time.sleep(0.05)
                with open(trace, "rb") as f:
                    while True:
                        chunk = f.read()
                        if chunk:
                            last = time.time()
                            raw.extend(chunk)
                            payload, used = itm_port0(raw)
                            del raw[:used]
                            if expiterations > 1 and start is not None:
                                pb.update(payload.count(b"+"))
                            text.extend(payload)
                            if start is None:
                                match = self.start_pat.match(text)
                                if match is not None:
                                    start = match.end()
                                    if expiterations > 1:
                                        pb.update(text[start:].count(b"+"))
                            if start is not None and b"#" in text[start:]:
                                break
                        elif proc.poll() is not None:
                            raise RuntimeError(f"OpenOCD exited with {proc.returncode}")
                        elif time.time() - last > self.timeout:
                            raise RuntimeError("Timeout waiting for ITM output")
                        else:
                            time.sleep(0.05)
            finally:
                proc.terminate()
                proc.wait()
        if expiterations > 1:
            pb.close()
        output = text[start:]
        return output[:output.index(b"#")].decode('utf-8', 'ignore')


class ChipWhisperer(mupq.Platform):

    # Start pattern is at least five equal signs
    start_pat = re.compile('.*={4,}\n', re.DOTALL)
    # End pattern is a hash with a newline
    end_pat = re.compile('.*#\n', re.DOTALL)

    def __init__(self):
        super().__init__()
        self.platformname = "cw"
        self.scope = cw.scope()
        self.target = cw.target(self.scope)
        self.scope.default_setup()

    def __enter__(self):
        return super().__enter__()

    def __exit__(self, *args, **kwargs):
        self.target.close()
        return super().__exit__(*args, **kwargs)

    def device(self):
        return self.wrapper

    def reset_target(self):
        self.scope.io.nrst = 'low'
        time.sleep(0.05)
        self.scope.io.nrst = 'high'
        time.sleep(0.05)

    def flash(self, binary_path):
        prog = cw.programmers.STM32FProgrammer()
        prog.scope = self.scope
        prog.open()
        prog.find()
        prog.erase()
        prog.program(binary_path, memtype="flash", verify=False)
        prog.close()

    def run(self, binary_path, expiterations=1):
        if expiterations > 1:
            pb = tqdm.tqdm(total=expiterations, leave=False, desc="Running...")
        self.flash(binary_path)
        self.target.flush()
        self.reset_target()
        data = ''
        # Wait for the first equal sign
        while '=' not in data:
            data += self.target.read()
        # Wait for the end of the equal delimiter
        match = None
        while match is None:
            buf = self.target.read()
            data += buf
            match = self.start_pat.match(data)
        # Remove the start pattern
        data = data[match.end():]
        # Wait for the end
        match = None
        while match is None:
            buf = self.target.read()
            data += buf
            if expiterations > 1:
                if "+" in buf:
                    pb.update(buf.count("+"))
                else:
                    pb.refresh()
            match = self.end_pat.match(data)
        # Remove stop pattern and return
        if expiterations > 1:
            pb.close()
        return data[:match.end() - 2]
