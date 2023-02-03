line_number_with_fields = 7
median_read_length_field = 17

# noinspection PyUnresolvedReferences
with open(snakemake.input[0], "r", encoding="utf-8") as f:
    for i, line in enumerate(f):
        if i == line_number_with_fields:
            # for some reason the median_read_length is not guaranteed to be
            # an integer, so we have to round it
            median_read_length = round(float(line.split()[median_read_length_field]))
            break

# noinspection PyUnresolvedReferences
with open(snakemake.output[0], "w", encoding="utf-8") as f:
    f.write(str(median_read_length))