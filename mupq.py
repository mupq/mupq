import abc
from collections import defaultdict
import re
import os
import logging
import subprocess
import hashlib
import time
from datetime import datetime


logging.basicConfig(level=logging.DEBUG)


class Implementation(object):
    """Contains some properties of a scheme implementation"""

    #: regex to parse the paths into schemes
    _path_regex = re.compile(
        r'(?P<type>crypto_sign|crypto_kem)/'
        r'(?P<scheme>\S+)/'
        r'(?P<implementation>\S+)/?$')

    def __init__(self, project, primitive, scheme, implementation, path):
        """Sets up this scheme"""
        self.project = project
        self.primitive = primitive
        self.scheme = scheme
        self.implementation = implementation
        self.path = path

    @classmethod
    def from_path(cls, project, path):
        """
        Construct a scheme implemenation from a path

        Specify the project that owns it
        """
        matches = cls._path_regex.match(path)
        if not matches:
            raise Exception(f"Unexpected path format: '{path}'")
        return cls(project,
                   matches.group("type"),
                   matches.group("scheme"),
                   matches.group("implementation"),
                   path)

    def get_binary_path(self, type_):
        return (f"bin/{self.project}_{self.primitive}_{self.scheme}"
                f"_{self.implementation}_{type_}.bin")

    def build_binary(self, type_):
        subprocess.check_call(
            ['make',
             f'PROJECT={self.project}',
             f"TYPE={self.primitive.split('_')[1]}",
             f"SCHEME={self.scheme}",
             f"IMPLEMENTATION={self.implementation}",
             self.get_binary_path(type_)])

    def __str__(self):
        return f"{self.scheme} - {self.implementation}"


class PlatformSettings(object):
    """Contains the settings for a certain platform"""
    scheme_folders = [
        ('pqclean', 'mupq/pqclean/crypto_kem'),
        ('pqclean', 'mupq/pqclean/crypto_sign'),
    ]
    skip_list = []
    name = None

    def __init__(self):
        self.log = logging.getLogger(__class__.__name__)

    def __str__(self):
        return self.name

    def get_implementations(self, all=False):
        """Get the schemes"""
        try:
            for (parent, scheme_folder) in self.scheme_folders:
                for scheme in os.listdir(scheme_folder):
                    scheme_path = os.path.join(scheme_folder, scheme)
                    if not os.path.isdir(scheme_path):
                        continue
                    for implementation_path in os.listdir(scheme_path):
                        path = os.path.join(scheme_path,
                                            implementation_path)
                        impl = Implementation.from_path(parent, path)
                        if not all and self.should_skip(impl):
                            continue
                        yield impl
        except FileNotFoundError as e:
            raise Exception(
                "There is no bin/ folder. Please first make binaries."
            ) from e

    def should_skip(self, impl):
        """Should this Implementation be skipped?"""
        for item in self.skip_list:
            match = len(item) > 0
            for attribute, value in item.items():
                if getattr(impl, attribute) != value:
                    match = False
            if match:
                return True
        return False


class Platform(abc.ABC):
    """Generic platform interface"""

    def __init__(self):
        self.log = logging.getLogger("platform interface")

    def device(self):
        raise NotImplementedError("Override this")

    @abc.abstractmethod
    def flash(self, binary_path):
        self.log.info("Flashing %s to device", binary_path)
        self.state = 'waiting'
        self.equals_seen = 0

    def run(self, binary_path):
        """Runs the flashed target and collects the result"""
        self.flash(binary_path)
        while not self._wait_for_start():
            pass
        self.log.info("Output started")
        return self._read_output()

    def _wait_for_start(self):
        """Waits until we read five equals signs"""
        while self.state == 'waiting':
            x = self.device().read()
            if x == b'':
                self.log.warning(
                    "timed out while waiting for the markers, reflashing")
                return False
            elif x == b'=':
                self.equals_seen += 1
                continue
            elif self.equals_seen > 5:
                self.state = 'beginning'
                self.log.debug("Found output marker")
            elif self.equals_seen > 1:
                self.log.warning(
                    "Got garbage after finding first equals sign, restarting"
                )
                return False
        # Read remaining = signs
        while self.state == 'beginning':
            x = self.device().read()
            # Consume remaining =
            if x != b'=':
                self.output = [x]
                self.state = 'reading'
                break
        return True

    def _read_output(self):
        while self.state == 'reading':
            x = self.device().read()
            if x == b'#':
                self.state = 'done'
                break
            elif x != b'':
                self.output.append(x)
        output = b''.join(self.output).decode('utf-8', 'ignore')
        # sometimes there's a line full of markers; strip out to avoid errors
        lines = (x for x in output.split('\n') if not all(c == '=' for c in x))
        return "{}\n".format('\n'.join(lines))


class BoardTestCase(abc.ABC):
    """
    Generic test class to run tests on all schemes.

    Generally you want to override run_test to parse the output of the tests
    running on the board.
    """
    test_type = 'undefined'

    def __init__(self, settings, interface):
        self.platform_settings = settings
        self.interface = interface
        self.log = logging.getLogger(__class__.__name__)

    def get_implementations(self, all=False):
        return self.platform_settings.get_implementations(all)

    @abc.abstractmethod
    def run_test(self, implementation):
        implementation.build_binary(f'{self.test_type}')
        binary = implementation.get_binary_path(f'{self.test_type}')
        return self.interface.run(binary)

    def test_all(self):
        for implementation in self.get_implementations():
            self.run_test(implementation)


class SimpleTest(BoardTestCase):
    test_type = 'test'

    def run_test(self, *args, **kwargs):
        output = super().run_test(*args, **kwargs).strip()
        if output.count("ERROR") or output.count("OK") != 30:
            self.log.error("Failed!")
        else:
            self.log.info("Success")


class StackBenchmark(BoardTestCase):
    test_type = 'stack'

    def run_test(self, implementation):
        self.log.info("Benchmarking %s", implementation)
        output = super().run_test(implementation)
        assert 'ERROR KEYS' not in output

        timestamp = datetime.fromtimestamp(
            time.time()).strftime('%Y%m%d%H%M%S')
        filename = os.path.join(
            'benchmarks/',
            self.test_type, implementation.primitive,
            implementation.scheme, implementation.implementation,
            timestamp)
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'w') as f:
            f.write(output.strip())
        self.log.info("Wrote benchmark!")


class SpeedBenchmark(StackBenchmark):
    test_type = 'speed'

class HashingBenchmark(StackBenchmark):
    test_type = 'hashing'

class TestVectors(BoardTestCase):
    test_type = 'testvectors'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.testvectorhash = dict()

    def hash_output(self, output):
        hash = hashlib.sha3_256(output.strip()).hexdigest()
        return hash

    def run_test(self, implementation):
        checksum = self.hash_output(
            super().run_test(implementation).encode('utf-8'))
        assert self.testvectorhash[implementation.scheme] == checksum

    def _prepare_testvectors(self):
        for scheme, implementations in self.schemes.items():
            for impl in implementations:
                if impl.implementation not in ('ref', 'clean'):
                    continue
                # Build host version
                self.log.info("Running %s on host", impl)
                binpath = impl.get_binary_path(self.test_type)
                hostbin = (binpath
                           .replace('bin/', 'bin-host/')
                           .replace('.bin', ''))
                subprocess.check_call(['make', hostbin])
                checksum = self.hash_output(
                        subprocess.check_output(
                            [hostbin],
                            stderr=subprocess.DEVNULL,
                        ))
                self.testvectorhash[scheme] = checksum
                break

    def test_all(self):
        self.schemes = defaultdict(list)
        for implementation in self.get_implementations(all=True):
            self.schemes[implementation.scheme].append(implementation)

        self._prepare_testvectors()

        for implementation in self.get_implementations():
            self.run_test(implementation)
