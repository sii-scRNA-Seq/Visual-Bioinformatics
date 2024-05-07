import { Injectable } from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class MockSocket {
  connect(): void { }

  emit(eventName: string, arg: any): void { }

  on(eventName: string, callback: any): void { }
}
