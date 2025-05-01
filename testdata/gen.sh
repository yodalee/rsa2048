# note that the testdata is irreproducible since we use ssh-keygen
ssh-keygen -t rsa -b 2048 -f out -m PEM
openssl rsa -in out -text -noout -inform PEM > rsa_key_openssl.txt
