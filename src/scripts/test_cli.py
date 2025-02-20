#!/usr/bin/python

"""
(C) 2018,2019 Jack Lloyd

Botan is released under the Simplified BSD License (see license.txt)
"""

import subprocess
import sys
import os
import logging
import optparse # pylint: disable=deprecated-module
import time
import shutil
import tempfile
import re
import random
import json
import binascii

# pylint: disable=global-statement,unused-argument

CLI_PATH = None
TESTS_RUN = 0
TESTS_FAILED = 0

class TestLogHandler(logging.StreamHandler, object):
    def emit(self, record):
        # Do the default stuff first
        super(TestLogHandler, self).emit(record)
        if record.levelno >= logging.ERROR:
            global TESTS_FAILED
            TESTS_FAILED += 1

def setup_logging(options):
    if options.verbose:
        log_level = logging.DEBUG
    elif options.quiet:
        log_level = logging.WARNING
    else:
        log_level = logging.INFO

    lh = TestLogHandler(sys.stdout)
    lh.setFormatter(logging.Formatter('%(levelname) 7s: %(message)s'))
    logging.getLogger().addHandler(lh)
    logging.getLogger().setLevel(log_level)

def random_port_number():
    return random.randint(1024, 65535)

def test_cli(cmd, cmd_options, expected_output=None, cmd_input=None, expected_stderr=None, use_drbg=True):
    global TESTS_RUN

    TESTS_RUN += 1

    opt_list = []

    if isinstance(cmd_options, str):
        opt_list = cmd_options.split(' ')
    elif isinstance(cmd_options, list):
        opt_list = cmd_options

    if use_drbg:
        fixed_drbg_seed = "802" * 32
        drbg_options = ['--rng-type=drbg', '--drbg-seed=' + fixed_drbg_seed]
    else:
        drbg_options = []

    cmdline = [CLI_PATH, cmd] + drbg_options + opt_list

    logging.debug("Executing '%s'" % (' '.join([CLI_PATH, cmd] + opt_list)))

    stdout = None
    stderr = None

    if cmd_input is None:
        proc = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
    else:
        proc = subprocess.Popen(cmdline, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate(cmd_input.encode())

    if stderr:
        if expected_stderr is None:
            logging.error("Got output on stderr %s (stdout was %s)", stderr, stdout)
        else:
            if stderr != expected_stderr:
                logging.error("Got output on stderr %s which did not match expected value %s", stderr, expected_stderr)

    output = stdout.decode('ascii').strip()

    if expected_output is not None:
        if output != expected_output:
            logging.error("Got unexpected output running cmd %s %s", cmd, cmd_options)
            logging.info("Output lengths %d vs expected %d", len(output), len(expected_output))
            logging.info("Got %s", output)
            logging.info("Exp %s", expected_output)

    return output

def check_for_command(cmd):
    cmdline = [CLI_PATH, 'has_command', cmd]
    proc = subprocess.Popen(cmdline)
    proc.communicate()

    return proc.returncode == 0

def cli_config_tests(_tmp_dir):
    prefix = test_cli("config", "prefix")
    cflags = test_cli("config", "cflags")
    ldflags = test_cli("config", "ldflags")
    libs = test_cli("config", "libs")

    if len(prefix) < 4 or prefix[0] != '/':
        logging.error("Bad prefix %s" % (prefix))
    if ("-I%s" % (prefix)) not in cflags:
        logging.error("Bad cflags %s" % (cflags))
    if ("-L%s" % (prefix)) not in ldflags:
        logging.error("Bad ldflags %s" % (ldflags))
    if "-lbotan-2" not in libs:
        logging.error("Bad libs %s" % (libs))

def cli_help_tests(_tmp_dir):
    output = test_cli("help", None, None)

    # Maybe test format somehow??
    if len(output) < 500:
        logging.error("Help output seems very short")

def cli_version_tests(_tmp_dir):
    output = test_cli("version", None, None)

    version_re = re.compile(r'[0-9]\.[0-9]+\.[0-9]')
    if not version_re.match(output):
        logging.error("Unexpected version output %s" % (output))

    output = test_cli("version", ["--full"], None, None)
    version_full_re = re.compile(r'Botan [0-9]\.[0-9]+\.[0-9] \(.* revision .*, distribution .*\)')
    if not version_full_re.match(output):
        logging.error("Unexpected version output %s" % (output))

def cli_is_prime_tests(_tmp_dir):
    test_cli("is_prime", "5", "5 is probably prime")
    test_cli("is_prime", "9", "9 is composite")
    test_cli("is_prime", "548950623407687320763", "548950623407687320763 is probably prime")

def cli_gen_prime_tests(_tmp_dir):
    test_cli("gen_prime", "64", "15568813029901363163")
    test_cli("gen_prime", "128", "287193909494025008847286845478788766073")

def cli_entropy_tests(_tmp_dir):
    output = test_cli("entropy", ["all"], None)

    status_re = re.compile('Polling [a-z0-9_]+ gathered [0-9]+ bytes in [0-9]+ outputs with estimated entropy [0-9]+')
    unavail_re = re.compile('Source [a-z0-9_]+ is unavailable')
    comp_re = re.compile('Sample from [a-z0-9_]+ was .* compressed from [0-9]+ bytes to [0-9]+ bytes')
    output_re = re.compile(r'[A-F0-9]+(...)?')

    status_next = True

    for line in output.split('\n'):
        if comp_re.match(line):
            continue

        if status_next:
            if status_re.match(line) is not None:
                status_next = False
            elif unavail_re.match(line) is not None:
                pass
            else:
                logging.error('Unexpected status line %s', line)
                status_next = False
        else:
            if output_re.match(line) is None:
                logging.error('Unexpected sample line %s', line)
            status_next = True

def cli_factor_tests(_tmp_dir):
    test_cli("factor", "97", "97: 97")
    test_cli("factor", "9753893489562389", "9753893489562389: 21433 455087644733")
    test_cli("factor", "12019502040659149507", "12019502040659149507: 3298628633 3643787579")

def cli_mod_inverse_tests(_tmp_dir):
    test_cli("mod_inverse", "97 802", "339")
    test_cli("mod_inverse", "98 802", "0")

def cli_base64_tests(_tmp_dir):
    test_cli("base64_enc", "-", "YmVlcyE=", "bees!")
    test_cli("base64_dec", "-", "bees!", "YmVlcyE=")

def cli_base32_tests(_tmp_dir):
    test_cli("base32_enc", "-", "MJSWK4ZB", "bees!")
    test_cli("base32_dec", "-", "bees!", "MJSWK4ZB")

def cli_base58_tests(_tmp_dir):
    test_cli("base58_enc", "-", "C6sRAr4", "bees!")
    test_cli("base58_dec", "-", "bees!", "C6sRAr4")

    test_cli("base58_enc", ["--check", "-"], "Cjv15cdjaBc", "F00F")
    test_cli("base58_dec", ["--check", "-"], "F00F", "Cjv15cdjaBc")

def cli_hex_tests(_tmp_dir):
    test_cli("hex_enc", "-", "6265657321", "bees!")
    test_cli("hex_dec", "-", "bees!", "6265657321")

def cli_hash_tests(_tmp_dir):
    test_cli("hash", "--algo=SHA-256",
             "E3B0C44298FC1C149AFBF4C8996FB92427AE41E4649B934CA495991B7852B855 -", "")

    test_cli("hash", "--algo=SHA-256",
             "BA7816BF8F01CFEA414140DE5DAE2223B00361A396177A9CB410FF61F20015AD -", "abc")

    test_cli("hash", ["--algo=SHA-256", "--format=base64"],
             "ungWv48Bz+pBQUDeXa4iI7ADYaOWF3qctBD/YfIAFa0= -", "abc")

    test_cli("hash", ["--algo=SHA-224", "--format=base58", "--no-fsname"],
             "MuGc8HkSVyJjfMjPM5UQikPToBTzNucEghcGLe", "abc")

    test_cli("hash", ["--algo=SHA-224", "--format=base58check", "--no-fsname"],
             "3MmfMqgrhemdVa9bDAGfooukbviWtKMBx2xauL2RsyAe", "abc")

def cli_hmac_tests(tmp_dir):
    key_file = os.path.join(tmp_dir, 'hmac.key')

    test_cli("rng", ["64", "--output=%s" % (key_file)], "")

    test_cli("hmac", ["--no-fsname", "--hash=SHA-384", key_file, key_file],
             "E3A8529377030B28A7DBDFC50DDEC8E4ECEFB6EA850D95EB785938CD3E3AFEF9EF8B08AF219C1496633193468AB755CB")

def cli_bcrypt_tests(_tmp_dir):
    test_cli("gen_bcrypt", "--work-factor=4 s3kr1t",
             "$2a$04$0.8G7o08XYwvBBWA3l0WUujtwoGZgGDzVSN8fNkNqXikcK4A3lHPS")

    test_cli("check_bcrypt", "s3kr1t $2a$04$gHX4Qg7pDSJuXiPXnmt8leyb.FFzX1Bv4rXwIj2cPSakJ8zNnhIka",
             "Password is valid")

    test_cli("check_bcrypt", "santa $2a$04$gHX4Qg7pDSJuXiPXnmt8leyb.FFzX1Bv4rXwIj2cPSakJ8zNnhIka",
             "Password is NOT valid")

def cli_argon2_tests(_tmp_dir):
    password = "s3kr1t"
    expected = "$argon2id$v=19$m=8,t=1,p=1$2A+I9q2+ZayxDDYC5n2YWw$/Lhx+Jbtlpw+Kxpskfv7+AKhBL/5ebalTJkVC1O5+1E"
    test_cli("gen_argon2", ['--mem=8', password], expected)
    test_cli("gen_argon2", ['--mem=8', '--t=1', password], expected)
    test_cli("gen_argon2", ['--mem=8', '--t=1', '--p=1', password], expected)

    test_cli("check_argon2", [password, expected], "Password is valid")
    test_cli("check_argon2", ["guessing", expected], "Password is NOT valid")

def cli_gen_dl_group_tests(_tmp_dir):

    pem = """-----BEGIN X9.42 DH PARAMETERS-----
MIIBJAKBgwTw7LQiLkXJsrgMVQxTPlWaQlYz/raZ+5RtIZe4YluQgRQGPFADLZ/t
TOYzuIzZJFOcdKtEtrVkxZRGSkjZwKFKLUD6fzSjoC2M2EHktK/y5HsvxBxL4tKr
q1ffbyPQi+iBLYTZAXygvxj2vWyrvA+/w4nbt1fStCHTDhWjLWqFpV9nAoGDAKzA
HUu/IRl7OiUtW/dz36gzEJnaYtz4ZtJl0FG8RJiOe02lD8myqW2sVzYqMvKD0LGx
x9fdSKC1G+aZ/NWtqrQjb66Daf7b0ddDx+bfWTWJ2dOtZd8IL2rmQQJm+JogDi9i
huVYFicDNQGzi+nEKAzrZ1L/VxtiSiw/qw0IyOuVtz8CFjgPiPatvmWssQw2AuZ9
mFvAZ/8wal0=
-----END X9.42 DH PARAMETERS-----"""

    test_cli("gen_dl_group", "--pbits=1043", pem)

    dsa_grp = """-----BEGIN X9.42 DH PARAMETERS-----
MIIBHgKBgQCyP1vosC/axliM2hmJ9EOSdd1zBkuzMP25CYD8PFkRVrPLr1ClSUtn
eXTIsHToJ7d7sRwtidQGW9BrvUEyiAWE06W/wnLPxB3/g2/l/P2EhbNmNHAO7rV7
ZVz/uKR4Xcvzxg9uk5MpT1VsxA8H6VEwzefNF1Rya92rqGgBTNT3/wKBgC7HLL8A
Gu3tqJxTk1iNgojjOiSreLn6ihA8R8kQnRXDTNtDKz996KHGInfMBurUI1zPM3xq
bHc0CvU1Nf87enhPIretzJcFgiCWrNFUIC25zPEjp0s3/ERHT4Bi1TABZ3j6YUEQ
fnnj+9XriKKHf2WtX0T4FXorvnKq30m934rzAhUAvwhWDK3yZEmphc7dwl4/J3Zp
+MU=
-----END X9.42 DH PARAMETERS-----"""

    test_cli("gen_dl_group", ["--type=dsa", "--pbits=1024"], dsa_grp)


def cli_key_tests(tmp_dir):

    pem = """-----BEGIN PRIVATE KEY-----
MIGEAgEAMBAGByqGSM49AgEGBSuBBAAKBG0wawIBAQQg2A+I9q2+ZayxDDYC5n2Y
W8Bn/zBm4D3mwS5qMwADRDehRANCAATwnDFqsjXL9SD/Rr1Vy4pb79PswXdQNZBN
mlLtJ5JvZ0/p6zP3x+Y9yPIrAR8L/acG5ItSrAKXzzuqQQZMv4aN
-----END PRIVATE KEY-----"""

    priv_key = os.path.join(tmp_dir, 'priv.pem')
    pub_key = os.path.join(tmp_dir, 'pub.pem')
    pub_der_key = os.path.join(tmp_dir, 'pub.der')
    enc_pem = os.path.join(tmp_dir, 'priv_enc.pem')
    enc_der = os.path.join(tmp_dir, 'priv_enc.der')
    ca_cert = os.path.join(tmp_dir, 'ca.crt')
    crt_req = os.path.join(tmp_dir, 'crt.req')
    user_cert = os.path.join(tmp_dir, 'user.crt')

    test_cli("keygen", ["--algo=ECDSA", "--params=secp256k1"], pem)

    test_cli("keygen", ["--algo=ECDSA", "--params=secp256r1", "--output=" + priv_key], "")

    test_cli("pkcs8", "--pub-out --output=%s %s" % (pub_key, priv_key), "")
    test_cli("pkcs8", "--pub-out --der-out --output=%s %s" % (pub_der_key, priv_key), "")

    test_cli("pkcs8", "--pass-out=foof --der-out --output=%s %s" % (enc_der, priv_key), "")
    test_cli("pkcs8", "--pass-out=foof --output=%s %s" % (enc_pem, priv_key), "")

    dec_pem = test_cli("pkcs8", ["--pass-in=foof", enc_pem], None)
    dec_der = test_cli("pkcs8", ["--pass-in=foof", enc_der], None)

    if dec_pem != dec_der:
        logging.error("Problem decrypting PKCS8 key")

    test_cli("fingerprint", ['--no-fsname', pub_key],
             "83:FC:67:87:30:C7:0C:9C:54:9A:E7:A1:FA:25:83:4C:77:A4:43:16:33:6D:47:3C:CE:4B:91:62:30:97:62:D4")

    test_cli("fingerprint", ['--no-fsname', pub_der_key],
             "83:FC:67:87:30:C7:0C:9C:54:9A:E7:A1:FA:25:83:4C:77:A4:43:16:33:6D:47:3C:CE:4B:91:62:30:97:62:D4")

    valid_sig = "nI4mI1ec14Y7nYUWs2edysAVvkob0TWpmGh5rrYWDA+/W9Fj0ZM21qJw8qa3/avAOIVBO6hoMEVmfJYXlS+ReA=="

    test_cli("sign", "--provider=base %s %s" % (priv_key, pub_key), valid_sig)

    test_cli("verify", [pub_key, pub_key, '-'],
             "Signature is valid", valid_sig)

    test_cli("verify", [pub_key, pub_key, '-'],
             "Signature is invalid",
             valid_sig.replace("G", "H"))

    test_cli("gen_self_signed",
             [priv_key, "CA", "--ca", "--country=VT",
              "--dns=ca.example", "--hash=SHA-384", "--output="+ca_cert],
             "")

    test_cli("cert_verify", ca_cert, "Certificate did not validate - Cannot establish trust")

    cert_info = test_cli("cert_info", ['--fingerprint', ca_cert], None)

    if cert_info.find('Subject: CN="CA",C="VT"') < 0:
        logging.error('Unexpected output for cert_info command %s', cert_info)
    if cert_info.find('Subject keyid: 69DD911C9EEE3400C67CBC3F3056CBE711BD56AF9495013F') < 0:
        logging.error('Unexpected output for cert_info command %s', cert_info)

    test_cli("gen_pkcs10", "%s User --output=%s" % (priv_key, crt_req))

    test_cli("sign_cert", "%s %s %s --output=%s" % (ca_cert, priv_key, crt_req, user_cert))

    test_cli("cert_verify", [user_cert, ca_cert],
             "Certificate passes validation checks")

    test_cli("cert_verify", user_cert,
             "Certificate did not validate - Certificate issuer not found")

def cli_xmss_sign_tests(tmp_dir):
    priv_key = os.path.join(tmp_dir, 'priv.pem')
    pub_key = os.path.join(tmp_dir, 'pub.pem')
    pub_key2 = os.path.join(tmp_dir, 'pub2.pem')
    msg = os.path.join(tmp_dir, 'input')
    sig1 = os.path.join(tmp_dir, 'sig1')
    sig2 = os.path.join(tmp_dir, 'sig2')

    test_cli("rng", ['--output=%s' % (msg)], "")
    test_cli("hash", ["--no-fsname", msg], "E3B0C44298FC1C149AFBF4C8996FB92427AE41E4649B934CA495991B7852B855")

    test_cli("keygen", ["--algo=XMSS", "--output=%s" % (priv_key)], "")
    test_cli("hash", ["--no-fsname", priv_key], "32397312E3FAC9D6396C55FEEFFF11EE195E2D2D5B34279D2544AF27763B0946")

    test_cli("pkcs8", "--pub-out --output=%s %s" % (pub_key, priv_key), "")
    test_cli("fingerprint", ['--no-fsname', pub_key],
             "E2:BE:C8:6D:CF:4B:5D:67:AB:A1:C1:F8:36:79:D5:3B:D8:17:D5:E3:5B:BE:29:08:03:7E:6E:07:27:4E:16:46")

    # verify the key is updated after each signature:
    test_cli("sign", [priv_key, msg, "--output=%s" % (sig1)], "")
    test_cli("verify", [pub_key, msg, sig1], "Signature is valid")
    test_cli("hash", ["--no-fsname", sig1], "04AF45451C7A9AF2D828E1AD6EC262E012436F4087C5DA6F32C689D781E597D0")
    test_cli("hash", ["--no-fsname", priv_key], "649E54D334F78A6AAAE34CFABF62121C74909D80E4DC2FA240A6EE1848526094")

    test_cli("sign", [priv_key, msg, "--output=%s" % (sig2)], "")
    test_cli("verify", [pub_key, msg, sig2], "Signature is valid")
    test_cli("hash", ["--no-fsname", sig2], "0785A6AD54CC7D01F2BE2BC6463A3EAA1159792E52210ED754992C5068E8F24F")
    test_cli("hash", ["--no-fsname", priv_key], "04483FA5367A7340F4BF6160FABD5742258009E05F9584E8D9732660B132608E")

    # private key updates, public key is unchanged:
    test_cli("pkcs8", "--pub-out --output=%s %s" % (pub_key2, priv_key), "")
    test_cli("fingerprint", ['--no-fsname', pub_key2],
             "E2:BE:C8:6D:CF:4B:5D:67:AB:A1:C1:F8:36:79:D5:3B:D8:17:D5:E3:5B:BE:29:08:03:7E:6E:07:27:4E:16:46")

def cli_pbkdf_tune_tests(_tmp_dir):
    if not check_for_command("pbkdf_tune"):
        return

    expected = re.compile(r'For (default|[1-9][0-9]*) ms selected Scrypt\([0-9]+,[0-9]+,[0-9]+\) using [0-9]+ MiB')

    output = test_cli("pbkdf_tune", ["--check", "1", "10", "50", "default"], None).split('\n')

    for line in output:
        if expected.match(line) is None:
            logging.error("Unexpected line '%s'" % (line))

    expected_pbkdf2 = re.compile(r'For (default|[1-9][0-9]*) ms selected PBKDF2\(HMAC\(SHA-256\),[0-9]+\)')

    output = test_cli("pbkdf_tune", ["--algo=PBKDF2(SHA-256)", "--check", "1", "10", "50", "default"], None).split('\n')

    for line in output:
        if expected_pbkdf2.match(line) is None:
            logging.error("Unexpected line '%s'" % (line))

    expected_argon2 = re.compile(r'For (default|[1-9][0-9]*) ms selected Argon2id\([0-9]+,[0-9]+,[0-9]+\)')

    output = test_cli("pbkdf_tune", ["--algo=Argon2id", "--check", "1", "10", "50", "default"], None).split('\n')

    for line in output:
        if expected_argon2.match(line) is None:
            logging.error("Unexpected line '%s'" % (line))

def cli_psk_db_tests(tmp_dir):
    if not check_for_command("psk_get"):
        return

    psk_db = os.path.join(tmp_dir, 'psk.db')
    db_key1 = "909"*32
    db_key2 = "451"*32

    test_cli("psk_set", [psk_db, db_key1, "name", "F00FEE"], "")
    test_cli("psk_set", [psk_db, db_key2, "name", "C00FEE11"], "")
    test_cli("psk_set", [psk_db, db_key1, "name2", "50051029"], "")

    test_cli("psk_get", [psk_db, db_key1, "name"], "F00FEE")
    test_cli("psk_get", [psk_db, db_key2, "name"], "C00FEE11")

    test_cli("psk_list", [psk_db, db_key1], "name\nname2")
    test_cli("psk_list", [psk_db, db_key2], "name")

def cli_compress_tests(tmp_dir):

    if not check_for_command("compress"):
        return

    input_file = os.path.join(tmp_dir, 'input.txt')
    output_file = os.path.join(tmp_dir, 'input.txt.gz')

    with open(input_file, 'w') as f:
        f.write("hi there")
        f.close()

    test_cli("compress", input_file)

    if not os.access(output_file, os.R_OK):
        logging.error("Compression did not created expected output file")

    is_py3 = sys.version_info[0] == 3

    output_hdr = open(output_file, 'rb').read(2)

    if is_py3:
        if output_hdr[0] != 0x1F or output_hdr[1] != 0x8B:
            logging.error("Did not see expected gzip header")
    else:
        if ord(output_hdr[0]) != 0x1F or ord(output_hdr[1]) != 0x8B:
            logging.error("Did not see expected gzip header")

    os.unlink(input_file)

    test_cli("decompress", output_file)

    if not os.access(input_file, os.R_OK):
        logging.error("Decompression did not created expected output file")

    recovered = open(input_file).read()
    if recovered != "hi there":
        logging.error("Decompression did not recover original input")

def cli_rng_tests(_tmp_dir):
    test_cli("rng", "10", "D80F88F6ADBE65ACB10C")
    test_cli("rng", "16", "D80F88F6ADBE65ACB10C3602E67D985B")
    test_cli("rng", "10 6", "D80F88F6ADBE65ACB10C\n1B119CC068AF")

    test_cli("rng", ['--format=base64', '10'], "2A+I9q2+ZayxDA==")
    test_cli("rng", ['--format=base58', '10'], "D93XRyVfxqs7oR")
    test_cli("rng", ['--format=base58check', '10'], "2NS1jYUq92TyGFVnhVLa")

    hex_10 = re.compile('[A-F0-9]{20}')

    for rng in ['system', 'auto', 'entropy']:
        output = test_cli("rng", ["10", '--%s' % (rng)], use_drbg=False)
        if output == "D80F88F6ADBE65ACB10C":
            logging.error('RNG produced DRBG output')
        if hex_10.match(output) is None:
            logging.error('Unexpected RNG output %s' % (output))

    has_rdrand = test_cli("cpuid", []).find(' rdrand ') > 0

    if has_rdrand:
        output = test_cli("rng", ["10", '--rdrand'], use_drbg=False)

        if output == "D80F88F6ADBE65ACB10C":
            logging.error('RDRAND produced DRBG output')
        if hex_10.match(output) is None:
            logging.error('Unexpected RNG output %s' % (output))


def cli_pk_workfactor_tests(_tmp_dir):
    test_cli("pk_workfactor", "1024", "80")
    test_cli("pk_workfactor", "2048", "111")
    test_cli("pk_workfactor", ["--type=rsa", "512"], "58")
    test_cli("pk_workfactor", ["--type=dl", "512"], "58")
    test_cli("pk_workfactor", ["--type=dl_exp", "512"], "128")

def cli_dl_group_info_tests(_tmp_dir):

    dl_output = re.compile('(P|G) = [A-F0-9]+')

    for bits in [1024, 1536, 2048, 3072, 4096, 6144, 8192]:
        output = test_cli("dl_group_info", "modp/ietf/%d" % (bits))
        lines = output.split('\n')

        if len(lines) != 2:
            logging.error('Unexpected output from dl_group_info')

        for l in lines:
            if not dl_output.match(l):
                logging.error('Unexpected output from dl_group_info')



def cli_ec_group_info_tests(_tmp_dir):

    # pylint: disable=line-too-long
    secp256r1_info = """P = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
A = FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC
B = 5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
N = FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
G = 6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5"""

    secp256r1_pem = """-----BEGIN EC PARAMETERS-----
MIHgAgEBMCwGByqGSM49AQECIQD/////AAAAAQAAAAAAAAAAAAAAAP//////////
/////zBEBCD/////AAAAAQAAAAAAAAAAAAAAAP///////////////AQgWsY12Ko6
k+ez671VdpiGvGUdBrDMU7D2O848PifSYEsEQQRrF9Hy4SxCR/i85uVjpEDydwN9
gS3rM6D0oTlF2JjClk/jQuL+Gn+bjufrSnwPnhYrzjNXazFezsu2QGg3v1H1AiEA
/////wAAAAD//////////7zm+q2nF56E87nKwvxjJVECAQE=
-----END EC PARAMETERS-----"""

    test_cli("ec_group_info", "secp256r1", secp256r1_info)
    test_cli("ec_group_info", "--pem secp256r1", secp256r1_pem)

def cli_cpuid_tests(_tmp_dir):
    cpuid_output = test_cli("cpuid", [])

    if not cpuid_output.startswith('CPUID flags: '):
        logging.error('Unexpected cpuid output %s' % (cpuid_output))

    flag_re = re.compile('[a-z0-9_]+')
    flags = cpuid_output[13:].split(' ')
    for flag in flags:
        if flag_re.match(flag) is None:
            logging.error('Unexpected CPUID flag name %s' % (flag))

def cli_cc_enc_tests(_tmp_dir):
    test_cli("cc_encrypt", ["8028028028028029", "pass"], "4308989841607208")
    test_cli("cc_decrypt", ["4308989841607208", "pass"], "8028028028028027")

def cli_cert_issuance_tests(tmp_dir):
    root_key = os.path.join(tmp_dir, 'root.key')
    root_crt = os.path.join(tmp_dir, 'root.crt')
    int_key = os.path.join(tmp_dir, 'int.key')
    int_crt = os.path.join(tmp_dir, 'int.crt')
    int_csr = os.path.join(tmp_dir, 'int.csr')
    leaf_key = os.path.join(tmp_dir, 'leaf.key')
    leaf_crt = os.path.join(tmp_dir, 'leaf.crt')
    leaf_csr = os.path.join(tmp_dir, 'leaf.csr')

    test_cli("keygen", ["--params=2048", "--output=" + root_key], "")
    test_cli("keygen", ["--params=2048", "--output=" + int_key], "")
    test_cli("keygen", ["--params=2048", "--output=" + leaf_key], "")

    test_cli("gen_self_signed",
             [root_key, "Root", "--ca", "--path-limit=2", "--output="+root_crt], "")

    test_cli("gen_pkcs10", "%s Intermediate --ca --output=%s" % (int_key, int_csr))
    test_cli("sign_cert", "%s %s %s --output=%s" % (root_crt, root_key, int_csr, int_crt))

    test_cli("gen_pkcs10", "%s Leaf --output=%s" % (leaf_key, leaf_csr))
    test_cli("sign_cert", "%s %s %s --output=%s" % (int_crt, int_key, leaf_csr, leaf_crt))

    test_cli("cert_verify" "%s %s %s" % (leaf_crt, int_crt, root_crt), "Certificate passes validation checks")

def cli_timing_test_tests(_tmp_dir):

    timing_tests = ["bleichenbacher", "manger",
                    "ecdsa", "ecc_mul", "inverse_mod", "pow_mod",
                    "lucky13sec3", "lucky13sec4sha1",
                    "lucky13sec4sha256", "lucky13sec4sha384"]

    output_re = re.compile('[0-9]+;[0-9];[0-9]+')

    for suite in timing_tests:
        output = test_cli("timing_test", [suite, "--measurement-runs=16", "--warmup-runs=3"], None).split('\n')

        for line in output:
            if output_re.match(line) is None:
                logging.error("Unexpected output in timing_test %s: %s", suite, line)

def cli_tls_ciphersuite_tests(_tmp_dir):
    policies = ['default', 'suiteb_128', 'suiteb_192', 'strict', 'all']

    versions = ['tls1.0', 'tls1.1', 'tls1.2']

    ciphersuite_re = re.compile('^[A-Z0-9_]+$')

    for policy in policies:
        for version in versions:

            if version != 'tls1.2' and policy != 'all':
                continue

            output = test_cli("tls_ciphers", ["--version=" + version, "--policy=" + policy], None).split('\n')

            for line in output:
                if ciphersuite_re.match(line) is None:
                    logging.error("Unexpected ciphersuite line %s", line)

def cli_asn1_tests(_tmp_dir):
    input_pem = """-----BEGIN BLOB-----
MCACAQUTBnN0cmluZzEGAQH/AgFjBAUAAAAAAAMEAP///w==
-----END BLOB------
"""

    expected = """d= 0, l=  32: SEQUENCE
  d= 1, l=   1:  INTEGER                                    05
  d= 1, l=   6:  PRINTABLE STRING                           string
  d= 1, l=   6:  SET
  d= 2, l=   1:   BOOLEAN                                   true
  d= 2, l=   1:   INTEGER                                   63
  d= 1, l=   5:  OCTET STRING                               0000000000
  d= 1, l=   4:  BIT STRING                                 FFFFFF"""

    test_cli("asn1print", "--pem -", expected, input_pem)

def cli_tls_socket_tests(tmp_dir):
    client_msg = b'Client message %d\n' % (random.randint(0, 2**128))
    server_port = random_port_number()

    priv_key = os.path.join(tmp_dir, 'priv.pem')
    ca_cert = os.path.join(tmp_dir, 'ca.crt')
    crt_req = os.path.join(tmp_dir, 'crt.req')
    server_cert = os.path.join(tmp_dir, 'server.crt')

    test_cli("keygen", ["--algo=ECDSA", "--params=secp256r1", "--output=" + priv_key], "")

    test_cli("gen_self_signed",
             [priv_key, "CA", "--ca", "--country=VT",
              "--dns=ca.example", "--hash=SHA-384", "--output="+ca_cert],
             "")

    test_cli("cert_verify", ca_cert, "Certificate did not validate - Cannot establish trust")

    test_cli("gen_pkcs10", "%s localhost --output=%s" % (priv_key, crt_req))

    test_cli("sign_cert", "%s %s %s --output=%s" % (ca_cert, priv_key, crt_req, server_cert))

    tls_server = subprocess.Popen([CLI_PATH, 'tls_server', '--max-clients=1',
                                   '--port=%d' % (server_port), server_cert, priv_key],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    time.sleep(.5)

    tls_client = subprocess.Popen([CLI_PATH, 'tls_client', 'localhost',
                                   '--port=%d' % (server_port), '--trusted-cas=%s' % (ca_cert)],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    time.sleep(.5)

    tls_client.stdin.write(client_msg)
    tls_client.stdin.flush()

    time.sleep(.5)

    (stdout, stderr) = tls_client.communicate()

    if stderr:
        logging.error("Got unexpected stderr output %s" % (stderr))

    if b'Handshake complete' not in stdout:
        logging.error('Failed to complete handshake: %s' % (stdout))

    if client_msg not in stdout:
        logging.error("Missing client message from stdout %s" % (stdout))

    tls_server.communicate()

def cli_tls_http_server_tests(tmp_dir):
    if not check_for_command("tls_http_server"):
        return

    try:
        from http.client import HTTPSConnection
    except ImportError:
        try:
            from httplib import HTTPSConnection
        except ImportError:
            return
    import ssl

    server_port = random_port_number()

    priv_key = os.path.join(tmp_dir, 'priv.pem')
    ca_cert = os.path.join(tmp_dir, 'ca.crt')
    crt_req = os.path.join(tmp_dir, 'crt.req')
    server_cert = os.path.join(tmp_dir, 'server.crt')

    test_cli("keygen", ["--algo=ECDSA", "--params=secp384r1", "--output=" + priv_key], "")

    test_cli("gen_self_signed",
             [priv_key, "CA", "--ca", "--country=VT",
              "--dns=ca.example", "--hash=SHA-384", "--output="+ca_cert],
             "")

    test_cli("gen_pkcs10", "%s localhost --output=%s" % (priv_key, crt_req))

    test_cli("sign_cert", "%s %s %s --output=%s" % (ca_cert, priv_key, crt_req, server_cert))

    tls_server = subprocess.Popen([CLI_PATH, 'tls_http_server', '--max-clients=2',
                                   '--port=%d' % (server_port), server_cert, priv_key],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    time.sleep(.5)

    context = ssl.create_default_context(cafile=ca_cert)
    conn = HTTPSConnection('localhost', port=server_port, context=context)
    conn.request("GET", "/")
    resp = conn.getresponse()

    if resp.status != 200:
        logging.error('Unexpected response status %d' % (resp.status))

    body = str(resp.read())

    if body.find('TLS negotiation with Botan 2.') < 0:
        logging.error('Unexpected response body')

    conn.request("POST", "/logout")
    resp = conn.getresponse()

    if resp.status != 405:
        logging.error('Unexpected response status %d' % (resp.status))

    if sys.version_info.major >= 3:
        rc = tls_server.wait(5) # pylint: disable=too-many-function-args
    else:
        rc = tls_server.wait()

    if rc != 0:
        logging.error("Unexpected return code from https_server %d", rc)

def cli_tls_proxy_tests(tmp_dir):
    # pylint: disable=too-many-locals,too-many-statements
    if not check_for_command("tls_proxy"):
        return

    try:
        from http.client import HTTPSConnection
    except ImportError:
        try:
            from httplib import HTTPSConnection
        except ImportError:
            return

    try:
        from http.server import HTTPServer, BaseHTTPRequestHandler
    except ImportError:
        try:
            from BaseHTTPServer import HTTPServer, BaseHTTPRequestHandler
        except ImportError:
            return

    import ssl
    import threading

    server_port = random_port_number()
    proxy_port = random_port_number()

    while server_port == proxy_port:
        proxy_port = random_port_number()

    priv_key = os.path.join(tmp_dir, 'priv.pem')
    ca_cert = os.path.join(tmp_dir, 'ca.crt')
    crt_req = os.path.join(tmp_dir, 'crt.req')
    server_cert = os.path.join(tmp_dir, 'server.crt')

    test_cli("keygen", ["--algo=ECDSA", "--params=secp384r1", "--output=" + priv_key], "")

    test_cli("gen_self_signed",
             [priv_key, "CA", "--ca", "--country=VT",
              "--dns=ca.example", "--hash=SHA-384", "--output="+ca_cert],
             "")

    test_cli("gen_pkcs10", "%s localhost --output=%s" % (priv_key, crt_req))

    test_cli("sign_cert", "%s %s %s --output=%s" % (ca_cert, priv_key, crt_req, server_cert))

    tls_proxy = subprocess.Popen([CLI_PATH, 'tls_proxy', str(proxy_port), '127.0.0.1', str(server_port),
                                  server_cert, priv_key, '--output=/tmp/proxy.err', '--max-clients=2'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    time.sleep(.5)

    server_response = binascii.hexlify(os.urandom(32))

    def run_http_server():
        class Handler(BaseHTTPRequestHandler):

            def do_GET(self): # pylint: disable=invalid-name
                self.send_response(200)
                self.end_headers()
                self.wfile.write(server_response)

        httpd = HTTPServer(('', server_port), Handler)
        httpd.serve_forever()

    http_thread = threading.Thread(target=run_http_server)
    http_thread.daemon = True
    http_thread.start()

    time.sleep(.5)

    context = ssl.create_default_context(cafile=ca_cert)

    for _i in range(2):
        conn = HTTPSConnection('localhost', port=proxy_port, context=context)
        conn.request("GET", "/")
        resp = conn.getresponse()

        if resp.status != 200:
            logging.error('Unexpected response status %d' % (resp.status))

        body = resp.read()

        if body != server_response:
            logging.error('Unexpected response from server %s' % (body))

    if sys.version_info.major >= 3:
        rc = tls_proxy.wait(5) # pylint: disable=too-many-function-args
    else:
        rc = tls_proxy.wait()

    if rc != 0:
        logging.error('Unexpected return code %d', rc)

def cli_trust_root_tests(tmp_dir):
    pem_file = os.path.join(tmp_dir, 'pems')
    dn_file = os.path.join(tmp_dir, 'dns')

    test_cli("trust_roots", ['--dn-only', '--output=%s' % (dn_file)], "")

    dn_re = re.compile('(.+=\".+\")(,.+=\".+\")')
    for line in open(dn_file):
        if dn_re.match(line) is None:
            logging.error("Unexpected DN line %s", line)

    test_cli("trust_roots", ['--output=%s' % (pem_file)], "")

def cli_tss_tests(tmp_dir):
    data_file = os.path.join(tmp_dir, 'data')

    exp_hash = "53B3C59276AE30EA7FD882268E80FD96AD80CC9FEB15F9FB940E7C4B5CF80B9E"

    test_cli("rng", ["32", "--output=%s" % (data_file)], "")
    test_cli("hash", ["--no-fsname", data_file], exp_hash)

    m = 3
    n = 5

    test_cli("tss_split", [str(m), str(n), data_file, "--share-prefix=%s/split" % (tmp_dir)], "")

    share_files = []

    for i in range(1, n+1):
        share = os.path.join(tmp_dir, "split%d.tss" % (i))
        if not os.access(share, os.R_OK):
            logging.error("Failed to create expected split file %s", share)
        share_files.append(share)

    rec5 = os.path.join(tmp_dir, "recovered_5")
    test_cli("tss_recover", share_files + ["--output=%s" % (rec5)], "")
    test_cli("hash", ["--no-fsname", rec5], exp_hash)

    rec4 = os.path.join(tmp_dir, "recovered_4")
    test_cli("tss_recover", share_files[1:] + ["--output=%s" % (rec4)], "")
    test_cli("hash", ["--no-fsname", rec4], exp_hash)

    rec3 = os.path.join(tmp_dir, "recovered_3")
    test_cli("tss_recover", share_files[2:] + ["--output=%s" % (rec3)], "")
    test_cli("hash", ["--no-fsname", rec3], exp_hash)

    rec2 = os.path.join(tmp_dir, "recovered_2")
    test_cli("tss_recover", share_files[3:] + ["--output=%s" % (rec2)], "", None,
             b'Error: Insufficient shares to do TSS reconstruction\n')


def cli_pk_encrypt_tests(tmp_dir):
    input_file = os.path.join(tmp_dir, 'input')
    ctext_file = os.path.join(tmp_dir, 'ctext')
    recovered_file = os.path.join(tmp_dir, 'recovered')
    rsa_priv_key = os.path.join(tmp_dir, 'rsa.priv')
    rsa_pub_key = os.path.join(tmp_dir, 'rsa.pub')

    test_cli("keygen", ["--algo=RSA", "--provider=base", "--params=2048", "--output=%s" % (rsa_priv_key)], "")

    key_hash = "891A3AA179639796B7A6348D2F1C3A8CC7E0FFED38BAE29143DF9B8A55391F28"
    test_cli("hash", ["--no-fsname", "--algo=SHA-256", rsa_priv_key], key_hash)

    test_cli("pkcs8", ["--pub-out", "%s/rsa.priv" % (tmp_dir), "--output=%s" % (rsa_pub_key)], "")

    # Generate a random input file
    test_cli("rng", ["10", "16", "32", "--output=%s" % (input_file)], "")

    # Because we used a fixed DRBG for each invocation the same ctext is generated each time
    rng_output_hash = "32F5E7B61357DE8397EFDA1E598379DFD5EE21767BDF4E2A435F05117B836AC6"
    ctext_hash = "5F45F360CF431C3E1BC126B1DB20CFE7A869AE7B67484A64F426A6349245EB51"

    test_cli("hash", ["--no-fsname", "--algo=SHA-256", input_file], rng_output_hash)

    # Encrypt and verify ciphertext is the expected value
    test_cli("pk_encrypt", [rsa_pub_key, input_file, "--output=%s" % (ctext_file)], "")
    test_cli("hash", ["--no-fsname", "--algo=SHA-256", ctext_file], ctext_hash)

    # Decrypt and verify plaintext is recovered
    test_cli("pk_decrypt", [rsa_priv_key, ctext_file, "--output=%s" % (recovered_file)], "")
    test_cli("hash", ["--no-fsname", "--algo=SHA-256", recovered_file], rng_output_hash)

def cli_uuid_tests(_tmp_dir):
    fixed_drbg_seed = "802" * 32

    test_cli("uuid", ['--rng-type=drbg', '--drbg-seed=' + fixed_drbg_seed], "D80F88F6-ADBE-45AC-B10C-3602E67D985B")

    uuid_re = re.compile(r'[0-9A-F]{8}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{12}')

    output = test_cli("uuid", [])

    if uuid_re.match(output) is None:
        logging.error('Bad uuid output %s' % (output))

def cli_tls_client_hello_tests(_tmp_dir):

    # pylint: disable=line-too-long
    chello = "16030100cf010000cb03035b3cf2457b864d7bef2a4b1f84fc3ced2b68d9551f3455ffdd305af277a91bb200003a16b816b716ba16b9cca9cca8c02cc030c02bc02fc0adc0acc024c00ac028c014c023c009c027c013ccaa009f009ec09fc09e006b003900670033010000680000000e000c000009676d61696c2e636f6d000500050100000000000a001a0018001d0017001a0018001b0019001c01000101010201030104000b00020100000d00140012080508040806050106010401050306030403001600000017000000230000ff01000100"

    output = test_cli("tls_client_hello", ["--hex", "-"], None, chello)

    output_hash = "8EBFC3205ACFA98461128FE5D081D19254237AF84F7DAF000A3C992C3CF6DE44"
    test_cli("hash", ["--no-fsname", "--algo=SHA-256", "-"], output_hash, output)

def cli_speed_tests(_tmp_dir):
    # pylint: disable=too-many-branches

    msec = 1

    output = test_cli("speed", ["--msec=%d" % (msec), "--buf-size=64,512", "AES-128"], None).split('\n')

    if len(output) % 4 != 0:
        logging.error("Unexpected number of lines for AES-128 speed test")

    # pylint: disable=line-too-long
    format_re = re.compile(r'^AES-128 .* buffer size [0-9]+ bytes: [0-9]+\.[0-9]+ MiB\/sec .*\([0-9]+\.[0-9]+ MiB in [0-9]+\.[0-9]+ ms\)')
    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    output = test_cli("speed", ["--msec=%d" % (msec), "ChaCha20", "SHA-256", "HMAC(SHA-256)"], None).split('\n')

    # pylint: disable=line-too-long
    format_re = re.compile(r'^.* buffer size [0-9]+ bytes: [0-9]+\.[0-9]+ MiB\/sec .*\([0-9]+\.[0-9]+ MiB in [0-9]+\.[0-9]+ ms\)')
    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    output = test_cli("speed", ["--msec=%d" % (msec), "AES-128/GCM"], None).split('\n')
    format_re_ks = re.compile(r'^AES-128/GCM\(16\) .* [0-9]+ key schedule/sec; [0-9]+\.[0-9]+ ms/op .*\([0-9]+ (op|ops) in [0-9]+ ms\)')
    format_re_cipher = re.compile(r'^AES-128/GCM\(16\) .* buffer size [0-9]+ bytes: [0-9]+\.[0-9]+ MiB\/sec .*\([0-9]+\.[0-9]+ MiB in [0-9]+\.[0-9]+ ms\)')
    for line in output:
        if format_re_ks.match(line) is None:
            if format_re_cipher.match(line) is None:
                logging.error('Unexpected line %s', line)

    pk_algos = ["ECDSA", "ECDH", "SM2", "ECKCDSA", "ECGDSA", "GOST-34.10",
                "DH", "DSA", "ElGamal", "Ed25519", "Curve25519", "NEWHOPE", "McEliece",
                "RSA", "XMSS"]

    output = test_cli("speed", ["--msec=%d" % (msec)] + pk_algos, None).split('\n')

    # ECDSA-secp256r1 106 keygen/sec; 9.35 ms/op 37489733 cycles/op (1 op in 9 ms)
    format_re = re.compile(r'^.* [0-9]+ ([A-Za-z ]+)/sec; [0-9]+\.[0-9]+ ms/op .*\([0-9]+ (op|ops) in [0-9]+ ms\)')
    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    # these all have a common output format
    math_ops = ['mp_mul', 'mp_div', 'mp_div10', 'modexp', 'random_prime', 'inverse_mod',
                'rfc3394', 'fpe_fe1', 'ecdsa_recovery', 'ecc_init', 'poly_dbl',
                'bn_redc', 'nistp_redc', 'ecc_mult', 'ecc_ops', 'os2ecp', 'primality_test',
                'bcrypt', 'passhash9', 'argon2']

    format_re = re.compile(r'^.* [0-9]+ /sec; [0-9]+\.[0-9]+ ms/op .*\([0-9]+ (op|ops) in [0-9]+(\.[0-9]+)? ms\)')
    for op in math_ops:
        output = test_cli("speed", ["--msec=%d" % (msec), op], None).split('\n')
        for line in output:
            if format_re.match(line) is None:
                logging.error("Unexpected line %s", line)

    output = test_cli("speed", ["--msec=%d" % (msec), "scrypt"], None).split('\n')

    format_re = re.compile(r'^scrypt-[0-9]+-[0-9]+-[0-9]+ \([0-9]+ MiB\) [0-9]+ /sec; [0-9]+\.[0-9]+ ms/op .*\([0-9]+ (op|ops) in [0-9]+ ms\)')

    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    output = test_cli("speed", ["--msec=%d" % (msec), "RNG"], None).split('\n')

    # ChaCha_RNG generate buffer size 1024 bytes: 954.431 MiB/sec 4.01 cycles/byte (477.22 MiB in 500.00 ms)
    format_re = re.compile(r'^.* generate buffer size [0-9]+ bytes: [0-9]+\.[0-9]+ MiB/sec .*\([0-9]+\.[0-9]+ MiB in [0-9]+\.[0-9]+ ms')
    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    # Entropy source rdseed output 128 bytes estimated entropy 0 in 0.02168 ms total samples 32
    output = test_cli("speed", ["--msec=%d" % (msec), "entropy"], None).split('\n')
    format_re = re.compile(r'^Entropy source [_a-z]+ output [0-9]+ bytes estimated entropy [0-9]+ in [0-9]+\.[0-9]+ ms .*total samples [0-9]+')
    for line in output:
        if format_re.match(line) is None:
            logging.error("Unexpected line %s", line)

    output = test_cli("speed", ["--msec=%d" % (msec), "--format=json", "AES-128"], None)

    json_blob = json.loads(output)
    if len(json_blob) < 2:
        logging.error("Unexpected size for JSON output")

    for b in json_blob:
        for field in ['algo', 'op', 'events', 'bps', 'buf_size', 'nanos']:
            if field not in b:
                logging.error('Missing field %s in JSON record %s' % (field, b))


def main(args=None):
    if args is None:
        args = sys.argv

    parser = optparse.OptionParser(
        formatter=optparse.IndentedHelpFormatter(max_help_position=50))

    parser.add_option('--verbose', action='store_true', default=False)
    parser.add_option('--quiet', action='store_true', default=False)

    (options, args) = parser.parse_args(args)

    setup_logging(options)

    if len(args) < 2:
        logging.error("Usage: ./cli_tests.py path_to_botan_cli [test_regex]")
        return 1

    if not os.access(args[1], os.X_OK):
        logging.error("Could not access/execute %s", args[1])
        return 2

    global CLI_PATH
    CLI_PATH = args[1]

    test_regex = None
    if len(args) == 3:
        try:
            test_regex = re.compile(args[2])
        except re.error as e:
            logging.error("Invalid regex: %s", str(e))
            return 1

    start_time = time.time()

    test_fns = [
        cli_argon2_tests,
        cli_asn1_tests,
        cli_base32_tests,
        cli_base58_tests,
        cli_base64_tests,
        cli_bcrypt_tests,
        cli_cc_enc_tests,
        cli_cert_issuance_tests,
        cli_compress_tests,
        cli_config_tests,
        cli_cpuid_tests,
        cli_dl_group_info_tests,
        cli_ec_group_info_tests,
        cli_entropy_tests,
        cli_factor_tests,
        cli_gen_dl_group_tests,
        cli_gen_prime_tests,
        cli_hash_tests,
        cli_help_tests,
        cli_hex_tests,
        cli_hmac_tests,
        cli_is_prime_tests,
        cli_key_tests,
        cli_xmss_sign_tests,
        cli_mod_inverse_tests,
        cli_pbkdf_tune_tests,
        cli_pk_encrypt_tests,
        cli_pk_workfactor_tests,
        cli_psk_db_tests,
        cli_rng_tests,
        cli_speed_tests,
        cli_timing_test_tests,
        cli_tls_ciphersuite_tests,
        cli_tls_client_hello_tests,
        cli_tls_http_server_tests,
        cli_tls_proxy_tests,
        cli_tls_socket_tests,
        cli_trust_root_tests,
        cli_tss_tests,
        cli_uuid_tests,
        cli_version_tests,
        ]

    for fn in test_fns:
        fn_name = fn.__name__

        if test_regex is not None:
            if test_regex.search(fn_name) is None:
                continue

        logging.info("Running %s" % (fn_name))

        start = time.time()
        tmp_dir = tempfile.mkdtemp(prefix='botan_cli_')
        try:
            fn(tmp_dir)
        except Exception as e: # pylint: disable=broad-except
            logging.error("Test %s threw exception: %s", fn_name, e)

        shutil.rmtree(tmp_dir)
        end = time.time()
        logging.debug("Ran %s in %.02f sec", fn_name, end-start)

    end_time = time.time()

    print("Ran %d tests with %d failures in %.02f seconds" % (
        TESTS_RUN, TESTS_FAILED, end_time - start_time))

    if TESTS_FAILED > 0:
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
