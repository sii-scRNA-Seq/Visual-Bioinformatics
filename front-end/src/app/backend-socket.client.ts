import { Inject, Injectable } from '@angular/core';
import { Observable, ReplaySubject } from 'rxjs';
import { Socket } from 'socket.io-client';

import { BackendSocketClientInterface } from './backend-socket.client.interface';
import { SOCKET } from './socket';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient implements BackendSocketClientInterface {

  private readonly response$: ReplaySubject<any> = new ReplaySubject<any> (1);
  readonly response: Observable<any> = this.response$.asObservable();

  constructor(@Inject(SOCKET) private socket: Socket) {
    this.socket.connect();
  }

  sendRequest(message: any): void {
    this.socket.emit('json', message);
  }

  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);
    });
  }
}