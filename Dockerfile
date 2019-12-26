FROM alpine

WORKDIR /usr/src/seqem
COPY seqem1.0/src/ .
RUN apk add --update alpine-sdk
RUN make
RUN mv seqem /usr/local/bin/seqem
