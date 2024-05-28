import { Inject, Injectable } from '@angular/core';
import { Socket } from 'socket.io-client';

import { Request } from './request';
import { SOCKET } from './socket';
import { UserIdService } from './user-id.service';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient {

  private userId: string | null = null;

  constructor(@Inject(SOCKET) private socket: Socket, private userIdService: UserIdService) {
    this.userIdService.userId.subscribe(
      (res) => { this.userId = res; },
    );
    this.socket.connect();
  }

  sendRequest(request: Request): void {
    this.socket.emit('json', request);
  }

  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      if (JSON.parse(msg).user_id == this.userId) {
        foo(msg);
      }
    });
    console.log(this.userId);
  }
}
