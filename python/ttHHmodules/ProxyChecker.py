# Attention: this works only for versions of Python below 3.7!

import os
import subprocess

class ProxyChecker(object):  # Explicitly inherit from object for Python 2 compatibility
    def __init__(self, path_proxy):
        self.path_proxy = path_proxy

    def check(self):
        # Check if the proxy file exists
        if not os.path.exists(self.path_proxy):
            error_msg = "  [ProxyChecker] [Error] : Proxy file does not exist. Please make the proxy with command [ voms-proxy-init --voms cms --valid 96:00 --out {0} ]".format(self.path_proxy)
            print(error_msg)
            raise IOError(error_msg)  # Using IOError instead of FileNotFoundError for Py2 compatibility

        # Check the validity time left of the proxy
        try:
            process = subprocess.Popen(
                ['voms-proxy-info', '-file', self.path_proxy, '--timeleft'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                error_msg = "  [ProxyChecker] Checking proxy validity: {0}".format(stderr.decode('utf-8').strip())
                print(error_msg)
                raise subprocess.CalledProcessError(
                    process.returncode, process.args, output=stdout, stderr=stderr
                )

            time_left = int(stdout.decode('utf-8').strip())
            if time_left <= 0:
                error_msg = "  [ProxyChecker] [Error] : Proxy has expired."
                print(error_msg)
                raise RuntimeError(error_msg)
            else:
                print("  [ProxyChecker] Proxy is valid for {0} seconds.".format(time_left))

        except Exception as e:
            print("  [ProxyChecker] [Error] Checking proxy validity: {0}".format(str(e)))
            raise
