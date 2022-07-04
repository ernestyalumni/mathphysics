def write_proto(filepath, proto_message):
    try:
        f = open(filepath, "wb")
        f.write(proto_message.SerializeToString())
        f.close()
    except IOError:
        print(filepath + ": Could not write file")


def read_proto(filepath, proto_message):
    try:
        f = open(filepath, "rb")
        proto_message.ParseFromString(f.read())
        f.close()
    except IOError:
        print(filepath + ": Could not open file. Creating a new one.")