filtering = open('/mnt/lustre/home/cusanovich/centipede/hg19_jack_centipede_sorted_pwms.bed','r')
filtered = open('/mnt/lustre/home/cusanovich/centipede/hg19_jack_centipede_sorted_pwms_clean.bed','w')

previous = []
for line in filtering:
    liner = line.strip().split()
    if previous == []:
        previous = liner
        continue
    if float(liner[6]) < 10:
        continue
    if liner[0:4] == previous[0:4]:
        previous[4] = str(max(float(previous[4]),float(liner[4])))
        previous[5] = "."
        continue
    print >> filtered, "\t".join(previous)
    previous = liner

filtering.close()
filtered.close()