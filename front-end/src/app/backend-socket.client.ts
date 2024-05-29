import { Inject, Injectable } from '@angular/core';
import { Socket } from 'socket.io-client';

import { Request } from './request';
import { SOCKET } from './socket';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient {

  constructor(@Inject(SOCKET) private socket: Socket) {
    this.socket.connect();
  }

  sendRequest(request: Request): void {
    this.socket.emit('json', request);
  }

  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);
    });
  }
}
