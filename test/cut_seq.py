def cut_seq(seq, cdsLoc, length):
    '''
    获取指定序列按照长度切割的结果
    '''
    seq_list = list(seq)
    seq = ''
    new_loc = cdsLoc
    count = len(seq_list)
    if count > length:
        if cdsLoc < length//2:
            seq = seq_list[0:length]
            new_loc = cdsLoc
        elif (count-cdsLoc) < length//2:
            seq = seq_list[0-length:]
            new_loc = length-(count-cdsLoc)
        else:
            left=cdsLoc-length//2-1
            seq = seq_list[left:left+length]
            new_loc = length//2+1
    else:
        seq = seq_list
        new_loc = cdsLoc
    return (''.join(seq), new_loc)

seq="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
sub_seq=cut_seq(seq,10,10)
print(sub_seq)