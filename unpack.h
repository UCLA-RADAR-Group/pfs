void unpack_pfs_2c2b (unsigned char *buf, char *outbuf, int bufsize);
void unpack_pfs_2c4b (unsigned char *buf, char *outbuf, int bufsize);
void unpack_pfs_2c8b (unsigned char *buf, char *outbuf, int bufsize);

void unpack_pfs_4c2b_rcp (unsigned char *buf, char *rcp, int bufsize);
void unpack_pfs_4c2b_lcp (unsigned char *buf, char *lcp, int bufsize);

void unpack_pfs_4c4b_rcp (unsigned char *buf, char *rcp, int bufsize);
void unpack_pfs_4c4b_lcp (unsigned char *buf, char *lcp, int bufsize);

void unpack_pfs_signed16bits(char *buf, float *outbuf, int bufsize);

