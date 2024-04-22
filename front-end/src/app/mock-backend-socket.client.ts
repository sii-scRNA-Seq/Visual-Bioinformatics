import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable, ReplaySubject } from 'rxjs';

import { Block } from './block.interface';
import { Output, RawOutput } from './output';
import { OutputServiceInterface } from './output.service.interface';
import { BackendSocketClientInterface } from './backend-socket.client.interface';

@Injectable({
  providedIn: 'root'
})
export class MockBackendSocketClient implements BackendSocketClientInterface {
  private readonly response$: ReplaySubject<any> = new ReplaySubject<any> (1);
  readonly response: Observable<any> = this.response$.asObservable();
  
  sendRequest(message: any): void {
    type RawOutput = { [key: string]: string };
    const res: RawOutput = {};
    if (message === "image") {
      res['img'] = 'Image text';
      res['alttext'] = 'Alt text';
    } else if (message === "end_connection") {
      res['end_connection'] = 'end_connection';
    } else if (message === "error") {
      res['error'] = 'error';
    } else if (message === "invalid response") {
      // do nothing
    } else {
      res['text'] = 'Dummy text';
    }
    this.response$.next(res);
  }
}