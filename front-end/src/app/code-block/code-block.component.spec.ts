import { ComponentFixture, TestBed } from '@angular/core/testing';
import { CodeBlockComponent } from './code-block.component';
import { By } from '@angular/platform-browser';
import { MatCardModule } from '@angular/material/card';
import { BlockService } from '../block.service';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CodeBlockComponent],
      imports: [MatCardModule],
      providers: [BlockService]
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
    component.blockId = 'LoadData'
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('RemoveBlock', () => {

    it ('should remove block when button is clicked', () => {
      
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, "removeBlock");
      const button = fixture.debugElement.query(By.css("button"));
      button.triggerEventHandler("click", {})
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith("LoadData");

    })
  })

});
