import { Injectable, isDevMode } from '@angular/core';
import { Observable, ReplaySubject } from 'rxjs';
import { Socket, io } from 'socket.io-client';

import { BackendSocketClientInterface } from './backend-socket.client.interface';

@Injectable({
  providedIn: 'root'
})
export class MockBackendSocketClient implements BackendSocketClientInterface {

  private socket: Socket;

  private readonly response$: ReplaySubject<any> = new ReplaySubject<any> (1);
  readonly response: Observable<any> = this.response$.asObservable();
  
  constructor() {
    const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
    // TODO: Are these timeouts required?
    this.socket = io(url, {
      autoConnect: false,
      timeout: 1000000,   // 1000 seconds,
      ackTimeout: 1000000 // 1000 seconds,
    });
    
    this.socket.connect();
  }
  
  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);  
    })
  }

  sendRequest(message: any): void {
    type RawOutput = { [key: string]: string };
    const res: RawOutput = {};
    if (message === 'image') {
      res['img'] = 'Image text';
      res['alttext'] = 'Alt text';
    } else if (message === 'end_connection') {
      res['end_connection'] = 'end_connection';
    } else if (message === 'error') {
      res['error'] = 'error';
    } else if (message === 'invalid response') {
      // do nothing
    } else {
      res['text'] = 'Dummy text';
    }
    this.response$.next(res);
  }
}